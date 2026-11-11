#pragma once
#include <atomic>
#include <chrono>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "actlab/concurrency.hpp"

namespace actlab {

// One logged sample. Fixed-size and trivially copyable so it can cross the SPSC queue.
struct LogRow {
  static constexpr int kCols = 16;
  double v[kCols]{};
};

// Streams LogRows to a CSV file from a background thread. The producer (simulation or
// control thread) only touches the lock-free queue. RAII: the destructor stops the writer,
// drains whatever is left and closes the file. Move-only because it owns a thread and a FILE*.
class CsvTelemetry {
 public:
  CsvTelemetry(const std::string& path, std::vector<std::string> header)
      : q_(std::make_unique<SpscQueue<LogRow, 1u << 15>>()) {
    file_ = std::fopen(path.c_str(), "w");
    if (!file_) throw std::runtime_error("cannot open " + path);
    ncols_ = static_cast<int>(header.size());
    if (ncols_ > LogRow::kCols) throw std::invalid_argument("too many columns");
    for (int i = 0; i < ncols_; ++i) std::fprintf(file_, "%s%s", header[i].c_str(), i + 1 < ncols_ ? "," : "\n");
    worker_ = std::thread([this] { run(); });
  }
  CsvTelemetry(const CsvTelemetry&) = delete;
  CsvTelemetry& operator=(const CsvTelemetry&) = delete;
  CsvTelemetry(CsvTelemetry&&) = delete;  // the worker captures `this`
  CsvTelemetry& operator=(CsvTelemetry&&) = delete;

  ~CsvTelemetry() {
    stop_.store(true, std::memory_order_release);
    if (worker_.joinable()) worker_.join();
    std::fclose(file_);
  }

  // Real-time producers use try_log and count drops; offline simulation uses log,
  // which yields until there is space so no sample is lost.
  bool try_log(const LogRow& r) {
    if (q_->try_push(r)) return true;
    dropped_.fetch_add(1, std::memory_order_relaxed);
    return false;
  }
  void log(const LogRow& r) {
    while (!q_->try_push(r)) std::this_thread::yield();
  }
  std::size_t dropped() const { return dropped_.load(std::memory_order_relaxed); }

 private:
  void run() {
    for (;;) {
      bool any = false;
      while (auto r = q_->try_pop()) {
        any = true;
        for (int i = 0; i < ncols_; ++i)
          std::fprintf(file_, "%.9g%s", r->v[i], i + 1 < ncols_ ? "," : "\n");
      }
      if (!any) {
        if (stop_.load(std::memory_order_acquire) && q_->empty()) return;
        std::this_thread::sleep_for(std::chrono::microseconds(200));
      }
    }
  }

  std::unique_ptr<SpscQueue<LogRow, 1u << 15>> q_;
  std::FILE* file_{nullptr};
  int ncols_{0};
  std::atomic<bool> stop_{false};
  std::atomic<std::size_t> dropped_{0};
  std::thread worker_;
};

}  // namespace actlab

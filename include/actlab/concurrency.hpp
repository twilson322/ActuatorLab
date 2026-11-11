#pragma once
#include <array>
#include <atomic>
#include <cstddef>
#include <optional>
#include <type_traits>

namespace actlab {

// Single-producer / single-consumer ring buffer. Wait-free on both sides, never allocates
// after construction, so the producer can be a real-time control thread.
// Capacity must be a power of two; one slot is left empty to tell full from empty.
template <typename T, std::size_t Capacity>
class SpscQueue {
  static_assert((Capacity & (Capacity - 1)) == 0, "Capacity must be a power of two");
  static_assert(std::is_trivially_copyable_v<T>, "keep payloads trivially copyable");

 public:
  bool try_push(const T& v) {
    const std::size_t head = head_.load(std::memory_order_relaxed);
    const std::size_t next = (head + 1) & (Capacity - 1);
    if (next == tail_.load(std::memory_order_acquire)) return false;  // full
    buf_[head] = v;
    head_.store(next, std::memory_order_release);  // publishes buf_[head]
    return true;
  }

  std::optional<T> try_pop() {
    const std::size_t tail = tail_.load(std::memory_order_relaxed);
    if (tail == head_.load(std::memory_order_acquire)) return std::nullopt;  // empty
    T v = buf_[tail];
    tail_.store((tail + 1) & (Capacity - 1), std::memory_order_release);
    return v;
  }

  bool empty() const {
    return tail_.load(std::memory_order_acquire) == head_.load(std::memory_order_acquire);
  }

 private:
  // Separate cache lines so producer and consumer do not false-share.
  alignas(64) std::atomic<std::size_t> head_{0};
  alignas(64) std::atomic<std::size_t> tail_{0};
  alignas(64) std::array<T, Capacity> buf_{};
};

// Lock-free triple buffer: "latest value" exchange between one writer and one reader.
// The writer never blocks, the reader always gets the newest complete sample, and buffer
// ownership moves through a single atomic index, so there is no torn read and no data race
// on the payload (ThreadSanitizer agrees, which it would not for a naive seqlock).
template <typename T>
class TripleBuffer {
 public:
  T& write_buffer() { return buf_[back_]; }
  void publish() {
    const unsigned prev = middle_.exchange(back_ | kFresh, std::memory_order_acq_rel);
    back_ = prev & kIndex;
  }
  // Returns true if a new sample arrived since the last call.
  bool update() {
    if (!(middle_.load(std::memory_order_relaxed) & kFresh)) return false;
    const unsigned prev = middle_.exchange(front_, std::memory_order_acq_rel);
    front_ = prev & kIndex;
    return true;
  }
  const T& read_buffer() const { return buf_[front_]; }

 private:
  static constexpr unsigned kFresh = 4u, kIndex = 3u;
  std::array<T, 3> buf_{};
  unsigned back_{0};                 // writer-owned
  unsigned front_{1};                // reader-owned
  std::atomic<unsigned> middle_{2};  // shared, low bits = index, bit 2 = fresh
};

}  // namespace actlab

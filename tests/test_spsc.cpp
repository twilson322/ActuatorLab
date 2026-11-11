#include <gtest/gtest.h>

#include <atomic>
#include <thread>

#include "actlab/concurrency.hpp"

using namespace actlab;

TEST(Concurrency, SpscPreservesOrderAcrossThreads) {
  static SpscQueue<std::uint64_t, 1024> q;
  constexpr std::uint64_t kN = 1'000'000;
  std::thread prod([] {
    for (std::uint64_t i = 1; i <= kN; ++i)
      while (!q.try_push(i)) std::this_thread::yield();
  });
  std::uint64_t expect = 1;
  while (expect <= kN) {
    if (auto v = q.try_pop()) {
      ASSERT_EQ(*v, expect);
      ++expect;
    }
  }
  prod.join();
  EXPECT_TRUE(q.empty());
}

TEST(Concurrency, TripleBufferNeverTearsAndIsMonotonic) {
  struct Msg { std::uint64_t seq; std::uint64_t check[7]; };
  static TripleBuffer<Msg> tb;
  std::atomic<bool> done{false};
  std::thread writer([&] {
    for (std::uint64_t i = 1; i <= 500000; ++i) {
      Msg& m = tb.write_buffer();
      m.seq = i;
      for (auto& c : m.check) c = i * 2654435761u;
      tb.publish();
    }
    done.store(true);
  });
  std::uint64_t last = 0, reads = 0;
  while (!done.load() || tb.update()) {
    if (tb.update()) {
      const Msg& m = tb.read_buffer();
      ASSERT_GT(m.seq, last);  // only ever newer samples
      for (auto c : m.check) ASSERT_EQ(c, m.seq * 2654435761u);  // never a mix of two writes
      last = m.seq;
      ++reads;
    }
  }
  writer.join();
  EXPECT_GT(reads, 0u);
}

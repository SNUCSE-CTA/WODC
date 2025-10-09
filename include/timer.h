/******************************************************************************
 * timer.h
 * record the elapsed time in microseconds (10^{-6} second)
 *****************************************************************************/

#ifndef TIMER_H_
#define TIMER_H_

#include <sys/time.h>

#include <cstdlib>

class Timer {
 public:
  Timer() { m_start = timestamp(); }
  void restart() { m_start = timestamp(); }
  uint64_t elapsed() {
    uint64_t time = timestamp() - m_start;
    return time;
  }

 private:
  uint64_t m_start;

  // Returns a timestamp ('now') in microseconds
  uint64_t timestamp() {
    struct timeval tp;
    gettimeofday(&tp, nullptr);
    return ((uint64_t)(tp.tv_sec)) * 1000000 + tp.tv_usec;
  }
};

#endif  // TIMER_H_

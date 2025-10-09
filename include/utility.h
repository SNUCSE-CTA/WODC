#ifndef UTILITY_H_
#define UTILITY_H_

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#define _BITSET_  // use bit set to represent the adjacency matrix
#define _STATISTIC_
#define _KERNEL_
#define _RECOLOR_

// #define NDEBUG
#include <cassert>

#ifdef _BITSET_
#define set_bit(array, pos) (((array)[(pos) >> 3]) |= (1 << ((pos) & 7)))
#define reverse_bit(array, pos) (((array)[(pos) >> 3]) ^= (1 << ((pos) & 7)))
#define test_bit1(array, pos) (((array)[(pos) >> 3]) & (1 << ((pos) & 7)))
#define test_bit(array, pos) ((((array)[(pos) >> 3]) >> ((pos) & 7)) & 1)
#else
#define set_bit(array, pos) (((array)[pos]) = 1)
#define reverse_bit(array, pos) (((array)[pos]) = 1 - ((array)[pos]))
#define test_bit(array, pos) ((array)[pos])
#endif

using ui = uint32_t;   // vertex type
using ept = uint32_t;  // edge pointer type; unsigned int can be used to
// process upto two billion undirected edges

constexpr ui INVALID = std::numeric_limits<ui>::max();

class Utility {
 public:
  static FILE *open_file(const char *file_name, const char *mode) {
    FILE *f = fopen(file_name, mode);
    if (f == nullptr) {
      printf("Can not open file: %s\n", file_name);
      exit(1);
    }

    return f;
  }

  static std::string integer_to_string(long long number) {
    std::vector<ui> sequence;
    if (number == 0) sequence.push_back(0);
    while (number > 0) {
      sequence.push_back(number % 1000);
      number /= 1000;
    }

    char buf[5];
    std::string res;
    for (ui i = sequence.size(); i > 0; i--) {
      if (i == sequence.size())
        sprintf(buf, "%u", sequence[i - 1]);
      else
        sprintf(buf, ",%03u", sequence[i - 1]);
      res += std::string(buf);
    }
    return res;
  }
};

#endif  // UTILITY_H_

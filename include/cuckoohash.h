#ifndef CUCKOOHASH_H_
#define CUCKOOHASH_H_

#include <immintrin.h>

#include "DAG.h"

/**
 * from the source code of "Maximal Defective Clique Enumeration"
 * (Dai et al., SIGMOD 2023)
 */

constexpr int32_t unfilled = -1;

namespace Cuckoo {

class CuckooHash {
 private:
  /* data */
  int32_t capacity;
  int32_t mask;
  int32_t size;
  int32_t buff_size = sizeof(int32_t);
  int32_t *hashtable;

  void rehash(int32_t **_table) {
    int32_t oldcapacity = capacity;
    mask = mask == 0 ? 1 : ((mask << 1) | 1);
    capacity = (mask + 1) * buff_size;
    int32_t *newhash = new int32_t[capacity];
    memset((newhash), unfilled, sizeof(int32_t) * capacity);
    for (int32_t i = 0; i < oldcapacity; ++i) {
      if ((*_table)[i] != unfilled) insert((*_table)[i], &newhash);
    }
    std::swap((*_table), newhash);
    delete[] newhash;
  }
  void insert(const int32_t &_u, int32_t **_table) {
    int32_t hs = hash1(_u);
    for (int32_t i = 0; i < buff_size; ++i) {
      if ((*_table)[hs * buff_size + i] == unfilled) {
        (*_table)[hs * buff_size + i] = _u;
        return;
      }
    }
    hs = hash2(_u);
    for (int32_t i = 0; i < buff_size; ++i) {
      if ((*_table)[hs * buff_size + i] == unfilled) {
        (*_table)[hs * buff_size + i] = _u;
        return;
      }
    }

    bool use_hash1 = true;
    int32_t u = _u;
    for (int32_t i = 0; i < mask; ++i) {
      int32_t replaced;
      if (use_hash1)
        hs = hash1(u);
      else
        hs = hash2(u);
      int32_t j = 0;
      for (; j < buff_size; ++j) {
        if ((*_table)[hs * buff_size + j] == unfilled) break;
      }
      if (buff_size == j) {
        replaced = std::move((*_table)[hs * buff_size]);
        j = 1;
        for (; j < buff_size; j++) {
          (*_table)[hs * buff_size + j - 1] =
              std::move((*_table)[hs * buff_size + j]);
        }
        (*_table)[hs * buff_size + j - 1] = u;
      } else {
        replaced = std::move((*_table)[hs * buff_size + j]);
        (*_table)[hs * buff_size + j] = u;
      }
      use_hash1 = hs == hash2(replaced);
      u = std::move(replaced);
      if (u == unfilled) return;
    }
    rehash(_table);
    insert(u, _table);
  }

  int32_t hash1(const int32_t &x) { return x & mask; }
  int32_t hash2(const int32_t &x) { return ~x & mask; }

 public:
  CuckooHash(/* args */) {
    capacity = 0;
    hashtable = NULL;
    mask = 0;
    size = 0;
  }
  ~CuckooHash() {
    if (hashtable) delete[] hashtable;
  }

  void reserve(int32_t _size) {
    if (capacity >= _size) return;
    mask = mask == 0 ? 1 : ((mask << 1) | 1);
    while (_size >= mask * buff_size) mask = (mask << 1) | 1;
    capacity = (mask + 1) * buff_size;
    if (hashtable) delete[] hashtable;
    hashtable = new int32_t[capacity];
    memset(hashtable, unfilled, sizeof(int32_t) * capacity);
  }

  void insert(const int32_t &_u) {
    if (find(_u)) return;
    insert(_u, &hashtable);
    size++;
  }

  bool find(const int32_t &_u) {
    int32_t hs1 = hash1(_u);
    int32_t hs2 = hash2(_u);

    assert(buff_size == 4 && sizeof(int32_t) == 4);
    __m128i cmp = _mm_set1_epi32(_u);
    __m128i b1 = _mm_load_si128((__m128i *)&hashtable[buff_size * hs1]);
    __m128i b2 = _mm_load_si128((__m128i *)&hashtable[buff_size * hs2]);
    __m128i flag =
        _mm_or_si128(_mm_cmpeq_epi32(cmp, b1), _mm_cmpeq_epi32(cmp, b2));

    return _mm_movemask_epi8(flag) != 0;
  }
  int32_t getcapacity() { return capacity; }
  int32_t getsize() { return size; }
  int32_t getmask() { return mask; }
  int32_t *gethashtable() { return hashtable; }
};

static std::vector<CuckooHash> cuhash;

static bool isNbr(int32_t u, int32_t v) {
  if (cuhash[u].getsize() <= 0)
    return false;
  else
    return cuhash[u].find(v);
}

static void initHash(const DAG *g) {
  if (cuhash.empty()) cuhash.resize(g->n);
  for (int u = 0; u < g->n; ++u) {
    int d = g->deg[0][u];
    cuhash[u].reserve(d);
    for (int j = 0; j < d; ++j) {
      int v = g->adj[g->cd[0][u] + j];
      cuhash[u].insert(v);
      if (u > v && !isNbr(v, u)) {
        printf("Can not find the nbr %d:%d\n", v, u);
        abort();
      }
    }
  }
}
}  // namespace Cuckoo

#endif  // CUCKOOHASH_H_

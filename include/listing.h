#ifndef LISTING_H_
#define LISTING_H_

#include <random>

#include "DAG.h"
#include "bucket_heap.h"
#include "cuckoohash.h"
#include "graph.h"

#ifdef NO_PRACTICAL
#define NO_COLOR_REDUCTION
#endif  // NO_PRACTICAL

using std::pair;
using std::vector;

static uint64_t num_recursive_calls = 0;
static uint64_t num_skip_leaf = 0;

static vector<ui> C;  // current clique
static DAG* g;
static BucketHeap* h;
static vector<ui> maximum_defective;
static ui k;
static ui mis = 0;

static ui max_depth;

bool inducing = true;

static uint64_t num_maximal = 0;

inline ui degidx(ui u) {
#ifdef LISTING_MAXIMAL
  return g->idx[u];
#else
#ifdef REDUCING
  return g->idx[u];
#else
  return 0;
#endif  // REDUCING
#endif  // LISTING_MAXIMAL
}

uint64_t num_type1 = 0;
uint64_t num_type2 = 0;
uint64_t num_pruned = 0;
uint64_t num_leaf = 0;

#ifdef PRINT_STAT
std::vector<uint64_t> vec_type1;
std::vector<uint64_t> vec_type2;
std::vector<uint64_t> vec_type3;
std::vector<uint64_t> vec_solution;
std::vector<uint64_t> vec_pruned;
inline void add_type1(ui l) {
  num_type1 += 1;
  vec_type1[l] += 1;
}
inline void add_type2(ui l) {
  num_type2 += 1;
  vec_type2[l] += 1;
}
inline void add_type3(ui l) { vec_type3[l] += 1; }
inline void add_solution(ui l) { vec_solution[l] += 1; }
inline void add_pruned(ui l) {
  num_pruned += 1;
  vec_pruned[l] += 1;
}
#else
inline void add_type1(ui l) { num_type1 += 1; }
inline void add_type2(ui l) { num_type2 += 1; }
inline void add_type3(ui l) {}
inline void add_solution(ui l) {}
inline void add_pruned(ui l) { num_pruned += 1; }
#endif  // PRINT_STAT

void expand_clique(ui u, bool is_maximal) {
  C.push_back(u);
  if (is_maximal) {
    if (C.size() > maximum_defective.size()) {
      num_maximal += 1;
#ifndef LISTING_MAXIMAL
      maximum_defective = C;
      std::cout << "Found a new " << k << "-defective clique of size "
                << C.size() << " (# missing edges: " << mis << ")" << std::endl;
#endif  // LISTING_MAXIMAL
      add_solution(C.size());
    }
  }
}

void shrink_clique() { C.pop_back(); }

/*
  set g->ssub[l], g->outdeg[l], g->color_cnt[l], g->outadj, g->color[l],
  g->num_color[l], g->num_distinct_color[l]
*/
void induce_and_coloring_dag1(ui l) {
  static std::vector<bool> color_vis(g->n);
  ui num_colors = 0;
  for (ui i = 0; i < h->cap(l); ++i) {
    g->ssub[l][i] = h->get(l, i);
    assert(!h->is_popped(g->ssub[l][i]));
  }
  std::sort(g->ssub[l], g->ssub[l] + h->cap(l));

  for (ui j = 0; j < h->cap(l); ++j) {
    ui v = g->ssub[l][j];
    assert(h->lab(v) == l && !h->is_popped(v));

    g->outdeg[l][v] = g->outdeg[l - 1][v];
    g->color_cnt[l][j] = 0;

    if (l == 1 || v != g->ssub[l - 1][j]) {
      for (ui s = g->outcd[v]; s < g->outcd[v] + g->outdeg[l][v]; ++s) {
        ui w = g->outadj[s];
        if (h->is_popped(w) || h->lab(w) < l) {
          std::swap(g->outadj[s--], g->outadj[g->outcd[v] + --g->outdeg[l][v]]);
        } else {
          color_vis[g->color[l][w]] = true;
        }
      }

      for (ui c = 0;; ++c) {
        if (color_vis[c] == false) {
          g->color[l][v] = c;
          num_colors = std::max(num_colors, c);
          g->color_cnt[l][c] += 1;
          break;
        }
      }

      for (ui s = g->outcd[v]; s < g->outcd[v] + g->outdeg[l][v]; ++s) {
        ui w = g->outadj[s];
        color_vis[g->color[l][w]] = false;
      }
    } else {
      g->color[l][v] = g->color[l - 1][v];
      num_colors = std::max(num_colors, g->color[l][v]);
      g->color_cnt[l][g->color[l][v]] += 1;
    }
  }
  if (h->cap(l) > 0) {
    num_colors += 1;
  }
  g->num_distinct_colors[l] = num_colors;
  g->num_colors[l] = num_colors;
}

void induce_and_coloring_dag2(ui l) {
  static std::vector<ui> color_vis1(g->n);
  static std::vector<ui> color_vis2(g->n);
  static std::vector<ui> anchor(g->n);
  static std::vector<std::vector<ui>> inadj(g->n);

  induce_and_coloring_dag1(l);

  ui prev_num_colors;
  do {
    if (g->num_colors[l] <= 3) {
      break;
    }
    prev_num_colors = g->num_colors[l];
    ui threshold = g->num_colors[l] - 2;
    ui num_colors = 0;

    for (ui j = 0; j < h->cap(l); ++j) {
      ui u = g->ssub[l][j];
      assert(h->lab(u) == l && !h->is_popped(u));

      inadj[u].resize(0);

      for (ui s = g->outcd[u]; s < g->outcd[u] + g->outdeg[l][u]; ++s) {
        ui v = g->outadj[s];
        assert(h->lab(v) == l && !h->is_popped(v));
        color_vis1[g->color[l][v]] += 1;
        anchor[g->color[l][v]] = v;
        inadj[v].push_back(u);
      }

      for (ui c = 0;; ++c) {
        if (color_vis1[c] == 0) {
          g->color[l][u] = c;
          break;
        }
      }

      if (g->color[l][u] > threshold) {
        bool updated = false;
        for (ui c = 0; c < threshold; ++c) {
          if (color_vis1[c] == 1) {
            const ui v = anchor[c];
            for (ui j = g->outcd[v]; j < g->outcd[v] + g->outdeg[l][v]; ++j) {
              const ui w = g->outadj[j];
              color_vis2[g->color[l][w]] = 1;
            }
            for (ui w : inadj[v]) {
              color_vis2[g->color[l][w]] = 1;
            }

            for (ui cc = 0; cc < threshold; ++cc) {
              if (c != cc && color_vis2[cc] == 0) {
                g->color[l][u] = c;
                g->color[l][v] = cc;
                color_vis1[c] = color_vis1[cc] = 0;
                color_vis2[c] = color_vis2[cc] = 0;
                num_colors = std::max(num_colors, cc);
                updated = true;
                break;
              }
            }

            for (ui j = g->outcd[v]; j < g->outcd[v] + g->outdeg[l][v]; ++j) {
              const ui w = g->outadj[j];
              color_vis2[g->color[l][w]] = 0;
            }
            for (ui w : inadj[v]) {
              color_vis2[g->color[l][w]] = 0;
            }

            if (updated) {
              break;
            }
          }
        }
      }

      num_colors = std::max(num_colors, g->color[l][u]);

      for (ui s = g->outcd[u]; s < g->outcd[u] + g->outdeg[l][u]; ++s) {
        ui v = g->outadj[s];
        color_vis1[g->color[l][v]] = 0;
      }
    }
    if (h->cap(l) > 0) {
      num_colors += 1;
    }
    g->num_distinct_colors[l] = num_colors;
    g->num_colors[l] = num_colors;
  } while (prev_num_colors < g->num_colors[l]);

  for (ui c = 0; c < g->num_colors[l]; ++c) {
    g->color_cnt[l][c] = 0;
  }
  for (ui j = 0; j < h->cap(l); ++j) {
    ui u = g->ssub[l][j];
    assert(g->color[l][u] < g->num_colors[l]);
    g->color_cnt[l][g->color[l][u]] += 1;
  }
}

void induce_and_coloring_dag(ui l) {
#ifdef TIGHT_RECOLORING
  induce_and_coloring_dag2(l);
#else
  induce_and_coloring_dag1(l);
#endif  // TIGHT_RECOLORING
}

void refine(ui u, const vector<pair<ui, ui>>& X, vector<pair<ui, ui>>& nX,
            vector<ui>& nX0, ui& xlab_max) {
#ifdef LISTING_MAXIMAL
  const ui l = C.size() + 1;
  nX.resize(0);
  nX0.resize(0);

  xlab_max = 0;

  for (auto [v, vl] : X) {
    if (l - vl <= k - mis + 1 && Cuckoo::isNbr(u, v)) {
      nX.emplace_back(v, vl + 1);
#ifndef NO_PIVOT_FROM_X
      if (vl == C.size()) {
        nX0.push_back(v);
      }
#endif  // NO_PIVOT_FROM_X
      xlab_max = std::max(xlab_max, vl + 1);
    } else if (l - vl <= k - mis) {
      nX.emplace_back(v, vl);
      xlab_max = std::max(xlab_max, vl);
    }
  }
#endif  // LISTING_MAXIMAL
}

static vector<ui> deg0;

void intersect1(ui u, ui num_cand, ui l, ui mu, ui* updated, ui& sz) {
  for (ui i = l + 1; i-- > 0 && l - i <= k - mis - mu;) {
    for (ui j = 0; j < h->cap(i); ++j) {
      const ui v = h->get(i, j);
      if (Cuckoo::isNbr(u, v) == (g->neg_idx[u] == 0)) {
        updated[sz++] = v;
      }
    }
#ifdef WORST_CASE_OPTIMAL
    if (i == l) {
      deg0[u] = sz;
    }
#endif  // WORST_CASE_OPTIMAL
  }
}

void intersect2(ui u, ui num_cand, ui l, ui mu, ui* updated, ui& sz) {
  for (ui j = g->cd[degidx(u)][u];
       j < g->cd[degidx(u)][u] + g->deg[degidx(u)][u]; ++j) {
    const ui v = g->adj[j];
    const ui vl = h->lab(v);
    if (mu + l - vl <= k - mis && !h->is_popped(v)) {
      updated[sz++] = v;
#ifdef WORST_CASE_OPTIMAL
      if (vl == l) {
        deg0[u] += 1;
      }
#endif  // WORST_CASE_OPTIMAL
    }
    assert(Cuckoo::isNbr(u, v) == (g->neg_idx[u] == 0));
  }
}

/* we are to add u to C and there are mu missing edges between u and C. */
void intersect(ui u, ui num_cand, ui l, ui mu, ui* updated, ui& sz) {
  sz = 0;
  deg0[u] = 0;
  if (num_cand * 1 < g->deg[degidx(u)][u]) {
    intersect1(u, num_cand, l, mu, updated, sz);
  } else {
    intersect2(u, num_cand, l, mu, updated, sz);
  }
  if (g->neg_idx[u] == 0) {
    if (num_cand < sz * 2) {
      ui cnt = num_cand - sz;
      sz = 0;
      deg0[u] = 0;
      if (cnt > 0) {
        for (ui i = l + 1; i-- > 0 && l - i <= k - mis - mu;) {
          for (ui j = 0; j < h->cap(i); ++j) {
            const ui v = h->get(i, j);
            if (u == v || !Cuckoo::isNbr(u, v)) {
              updated[sz++] = v;
              if (sz == cnt) {
#ifdef WORST_CASE_OPTIMAL
                if (i == l) {
                  deg0[u] = sz;
                }
#endif  // WORST_CASE_OPTIMAL
                goto OUT;
              }
            }
          }
#ifdef WORST_CASE_OPTIMAL
          if (i == l) {
            deg0[u] = sz;
          }
#endif  // WORST_CASE_OPTIMAL
        }
      }
    OUT:
      g->neg_idx[u] += 1;
      assert(cnt == sz);
    }
  } else {
    g->neg_idx[u] += 1;
  }
}

/**
 * Return 1 if pivot with no branch found,
 * Return 0 otherwise.
 */

vector<ui> num_removed;

vector<ui> prank;

std::mt19937 gen(2);

ui induce(ui num_cand, const vector<ui>& nX0, ui& pivot) {
  pivot = INVALID;
  if (inducing) {
    const ui l = C.size() + 1;
    num_removed.push_back(0);
    ui sz;
    for (auto v : nX0) {
      if (g->adj_off + std::min(num_cand, g->deg[degidx(v)][v]) >=
          g->adj_size) {
        g->resize_adj();
      }
      intersect(v, num_cand, l, 0, g->adj + g->adj_off, sz);
      g->idx[v] += 1;
      g->deg[degidx(v)][v] = sz;
      g->cd[degidx(v)][v] = g->adj_off;
      g->adj_off += sz;
      num_removed.back() += sz;
      if (g->neg_idx[v] == 0) {
        prank[v] = sz;
      } else {
        prank[v] = num_cand - sz;
        deg0[v] = h->cap(l) - deg0[v];
      }
#ifndef NDEBUG
      ui d = 0;
      for (ui i = 0; i < h->cap(l); ++i) {
        ui w = h->get(l, i);
        if (Cuckoo::isNbr(v, w)) {
          d += 1;
        }
      }
      assert(d == deg0[v]);
#endif  // NDEBUG
#ifdef WORST_CASE_OPTIMAL
      // if (pivot == INVALID || deg0[v] > deg0[pivot] ||
      //     (deg0[v] == deg0[pivot] && prank[v] > prank[pivot])) {
      //   pivot = v;
      // }
      if (deg0[v] == h->cap(l) &&
          (pivot == INVALID || prank[pivot] < prank[v])) {
        pivot = v;
      }
      // if (prank[v] == num_cand) {
      //   pivot = v;
      // }
#else
      if (pivot == INVALID || prank[v] > prank[pivot]) {
        pivot = v;
      }
#endif  // WORST_CASE_OPTIMAL
      if (prank[v] == num_cand) {
        return 1;
      }
    }

    for (ui j = 0; j < h->cap(l); ++j) {
      ui v = h->get(l, j);
      if (g->adj_off + std::min(num_cand, g->deg[degidx(v)][v]) >=
          g->adj_size) {
        g->resize_adj();
      }
      intersect(v, num_cand, l, 0, g->adj + g->adj_off, sz);
      g->idx[v] += 1;
      g->deg[degidx(v)][v] = sz;
      g->cd[degidx(v)][v] = g->adj_off;
      g->adj_off += sz;
      num_removed.back() += sz;
      if (g->neg_idx[v] == 0) {
        prank[v] = sz;
      } else {
        prank[v] = num_cand - sz;
        deg0[v] = h->cap(l) - deg0[v];
      }
#ifndef NDEBUG
      ui d = 0;
      for (ui i = 0; i < h->cap(l); ++i) {
        ui w = h->get(l, i);
        if (Cuckoo::isNbr(v, w)) {
          d += 1;
        }
      }
      assert(d == deg0[v]);
#endif  // NDEBUG
#ifndef RANDOM_PIVOT
#ifdef WORST_CASE_OPTIMAL
      if (pivot == INVALID || deg0[v] > deg0[pivot] ||
          (deg0[v] == deg0[pivot] && prank[v] > prank[pivot])) {
        pivot = v;
      }
#else
      if (pivot == INVALID || prank[v] > prank[pivot]) {
        pivot = v;
      }
#endif  // WORST_CASE_OPTIMAL
#endif  // RANDOM_PIVOT
    }
#ifdef RANDOM_PIVOT
    if (h->cap(l) > 0) {
      std::uniform_int_distribution<int> dis(0, h->cap(l) - 1);
      pivot = h->get(l, dis(gen));
    }
#endif  // RANDOM_PIVOT
#ifndef PIVOTING
    pivot = INVALID;
#endif  // PIVOTING
  }
  return 0;
}

void uninduce(ui num_cand, const vector<ui>& nX0) {
  if (inducing) {
    const ui l = C.size() + 1;
    ui n = num_removed.back();
    num_removed.pop_back();
    g->adj_off -= n;

    for (auto v : nX0) {
      ui prnk;
      if (g->neg_idx[v] == 0) {
        prnk = g->deg[degidx(v)][v];
      } else {
        prnk = num_cand - g->deg[degidx(v)][v];
      }
      g->idx[v] -= 1;
      if (g->neg_idx[v] > 0) {
        g->neg_idx[v] -= 1;
      }
      if (prnk == num_cand) {
        return;
      }
    }

    for (ui j = 0; j < h->cap(l); ++j) {
      ui v = h->get(l, j);
      g->idx[v] -= 1;
      if (g->neg_idx[v] > 0) {
        g->neg_idx[v] -= 1;
      }
    }
  }
}

ui max_mis(ui l) {
  ui m = h->cap(l);
  ui cur_mis = mis;

  for (ui i = l; i-- > 0 && cur_mis + l - i <= k;) {
    ui c = l - i;
    ui a = std::min(h->cap(i), (k - cur_mis) / c);
    m += a;
    cur_mis += a * c;
    assert(cur_mis <= k);
  }
  return m;
}

namespace club {
vector<vector<ui>> color_cnt;
vector<ui> total_color_cnt;
vector<ui> mis_cnt;
vector<ui> num_colors;
vector<vector<ui>> colors;
};  // namespace club

ui max_mis_defec(ui l, ui start) {
  ui m = 0;
  ui cur_mis = mis;

  for (ui i = start + 1; i-- > 0 && cur_mis + l - i <= k;) {
    ui c = l - i;
    ui a = std::min(h->cap(i), (k - cur_mis) / c);
    m += a;
    cur_mis += a * c;
    assert(cur_mis <= k);
  }

  if (l + m <= maximum_defective.size()) {
    return m;
  } else {
    // coloring-based upper bound
    m = 0;
    cur_mis = mis;
    static vector<ui> vis_colors;

    for (ui a = l - start; a <= k - cur_mis; ++a) {
      ui max_added = (k - cur_mis) / a;
      if (l + m + max_added <= maximum_defective.size()) {
        // pruned
        goto EXIT;
      }
      m += std::min(max_added, club::mis_cnt[a]);
      cur_mis += std::min(max_added, club::mis_cnt[a]) * a;
      assert(cur_mis <= k);

      if (cur_mis + a > k) {
        goto EXIT;
      }

      if (l == a) {
        for (ui c = 0; c < club::num_colors[0]; ++c) {
          for (ui j = 0; a + j <= k - cur_mis; ++j) {
            if (club::total_color_cnt[c] == 0) {
              m += 1;
              cur_mis += a;
              assert(cur_mis <= k);

              if (cur_mis + a > k) {
                goto EXIT;
              }
              vis_colors.push_back(c);
            } else if (a + club::total_color_cnt[c] <= k - cur_mis) {
              club::mis_cnt[a + club::total_color_cnt[c]] += 1;
            } else {
              break;
            }
            club::total_color_cnt[c] += 1;
          }
        }

        for (ui aa = a + 1; aa <= k - cur_mis; ++aa) {
          ui max_added = (k - cur_mis) / aa;
          if (l + m + max_added <= maximum_defective.size()) {
            // pruned
            goto EXIT;
          }
          m += std::min(max_added, club::mis_cnt[aa]);
          cur_mis += std::min(max_added, club::mis_cnt[aa]) * aa;
          assert(cur_mis <= k);

          if (cur_mis + aa > k) {
            goto EXIT;
          }
        }
        goto EXIT;
      } else {
        for (ui j = 0; j < h->cap(l - a); ++j) {
          ui u = h->get(l - a, j);
          ui uc = g->get_static_color(u);
          if (club::total_color_cnt[uc] == 0) {
            m += 1;
            cur_mis += a;
            assert(cur_mis <= k);

            if (cur_mis + a > k) {
              goto EXIT;
            }
            vis_colors.push_back(uc);
          } else if (a + club::total_color_cnt[uc] <= k - cur_mis) {
            club::mis_cnt[a + club::total_color_cnt[uc]] += 1;
          }
          club::total_color_cnt[uc] += 1;
        }
      }
    }
  EXIT:
    for (ui j = 2; j <= k - mis; ++j) {
      club::mis_cnt[j] = 0;
    }
    for (ui c : vis_colors) {
      club::total_color_cnt[c] = 0;
    }
    vis_colors.resize(0);
    assert(club::mis_cnt[1] == 0);

    return m;
  }
}

ui max_mis_color_dag(ui l, bool left_branch) {
#ifndef RECOLOR
  return g->n;
#else
  if (left_branch) {
    induce_and_coloring_dag(l);
  }
#ifndef NDEBUG
  ui ncolors0 = 0;
  static vector<bool> vis(g->n);
  static vector<ui> cnt0(g->n);

  for (ui j = 0; j < h->cap(l); ++j) {
    ui u = h->get(l, j);
    ui uc = g->color[l][u];
    if (vis[uc] == false) {
      vis[uc] = true;
      ncolors0 += 1;
    }
    cnt0[uc] += 1;
  }
  assert(ncolors0 == g->num_distinct_colors[l]);
  for (ui c = 0; c < g->num_colors[l]; ++c) {
    assert(g->color_cnt[l][c] == cnt0[c]);
  }
  for (ui j = 0; j < h->cap(l); ++j) {
    ui u = h->get(l, j);
    ui uc = g->color[l][u];
    vis[uc] = false;
    cnt0[uc] = 0;
  }
#endif  // NDEBUG

  ui m = g->num_distinct_colors[l];
  ui cur_mis = mis;

  if (l + m > maximum_defective.size() ||
      l + m + k - mis <= maximum_defective.size()) {
    return m + k - mis;
  } else {
    static vector<ui> vis_colors;

    for (ui c = 0; c < g->num_colors[l]; ++c) {
      if (g->color_cnt[l][c] > 1) {
        m += 1;
        cur_mis += 1;
        if (l + m > maximum_defective.size() || cur_mis == k) {
          goto EXIT;
        }
        for (ui cnt = 2; cnt <= std::min(g->color_cnt[l][c] - 1, k - cur_mis);
             ++cnt) {
          club::mis_cnt[cnt] += 1;
        }
      }
    }

    for (ui a = 1; a <= k - cur_mis; ++a) {
      ui max_added = (k - cur_mis) / a;
      if (l + m + max_added <= maximum_defective.size()) {
        // pruned
        goto EXIT;
      }
      m += std::min(max_added, club::mis_cnt[a]);
      cur_mis += std::min(max_added, club::mis_cnt[a]) * a;
      assert(cur_mis <= k);

      if (cur_mis + a > k) {
        goto EXIT;
      }

      if (l == a) {
        for (ui c = 0; c < club::num_colors[0]; ++c) {
          for (ui j = 0; a + j <= k - cur_mis; ++j) {
            if (club::total_color_cnt[c] == 0) {
              m += 1;
              cur_mis += a;
              assert(cur_mis <= k);

              if (cur_mis + a > k) {
                goto EXIT;
              }
              vis_colors.push_back(c);
            } else if (a + club::total_color_cnt[c] <= k - cur_mis) {
              club::mis_cnt[a + club::total_color_cnt[c]] += 1;
            } else {
              break;
            }
            club::total_color_cnt[c] += 1;
          }
        }

        for (ui aa = a + 1; aa <= k - cur_mis; ++aa) {
          ui max_added = (k - cur_mis) / aa;
          if (l + m + max_added <= maximum_defective.size()) {
            // pruned
            goto EXIT;
          }
          m += std::min(max_added, club::mis_cnt[aa]);
          cur_mis += std::min(max_added, club::mis_cnt[aa]) * aa;
          assert(cur_mis <= k);

          if (cur_mis + aa > k) {
            goto EXIT;
          }
        }
        goto EXIT;
      } else {
        for (ui j = 0; j < h->cap(l - a); ++j) {
          ui u = h->get(l - a, j);
          ui uc = g->get_static_color(u);
          if (club::total_color_cnt[uc] == 0) {
            m += 1;
            cur_mis += a;
            assert(cur_mis <= k);

            if (cur_mis + a > k) {
              goto EXIT;
            }
            vis_colors.push_back(uc);
          } else if (a + club::total_color_cnt[uc] <= k - cur_mis) {
            club::mis_cnt[a + club::total_color_cnt[uc]] += 1;
          }
          club::total_color_cnt[uc] += 1;
        }
      }
    }
  EXIT:
    for (ui j = 2; j <= k - mis; ++j) {
      club::mis_cnt[j] = 0;
    }
    for (ui c : vis_colors) {
      club::total_color_cnt[c] = 0;
    }
    vis_colors.resize(0);
    assert(club::mis_cnt[1] == 0);

    return m;
  }
#endif  // RECOLOR
}

ui max_mis_color(ui l, bool left_branch) {
  if (left_branch) {
    for (ui c : club::colors[l]) {
      club::color_cnt[l][c] = 0;
    }
    club::colors[l].resize(0);
    club::num_colors[l] = 0;
    for (ui j = 0; j < h->cap(l); ++j) {
      ui u = h->get(l, j);
      ui uc = g->get_static_color(u);
      if (club::color_cnt[l][uc] == 0) {
        club::num_colors[l] += 1;
        club::colors[l].push_back(uc);
      }
      club::color_cnt[l][uc] += 1;
    }
  }
#ifndef NDEBUG
  ui ncolors0 = 0;
  static vector<bool> vis(g->n);
  for (ui j = 0; j < h->cap(l); ++j) {
    ui u = h->get(l, j);
    ui uc = g->get_static_color(u);
    if (vis[uc] == false) {
      vis[uc] = true;
      ncolors0 += 1;
    }
  }
  for (ui j = 0; j < h->cap(l); ++j) {
    ui u = h->get(l, j);
    ui uc = g->get_static_color(u);
    vis[uc] = false;
  }
  assert(ncolors0 == club::num_colors[l]);
#endif  // NDEBUG

  ui m = club::num_colors[l];
  ui cur_mis = mis;

  if (l + m > maximum_defective.size() ||
      l + m + k - mis <= maximum_defective.size()) {
    return m + k - mis;
  } else {
    static vector<ui> vis_colors;

    for (ui c = 0; c < club::num_colors[0]; ++c) {
      if (club::color_cnt[l][c] > 1) {
        m += 1;
        cur_mis += 1;
        if (l + m > maximum_defective.size() || cur_mis == k) {
          goto EXIT;
        }
        for (ui cnt = 2;
             cnt <= std::min(club::color_cnt[l][c] - 1, k - cur_mis); ++cnt) {
          club::mis_cnt[cnt] += 1;
        }
      }
      if (club::color_cnt[l][c] > 0) {
        club::total_color_cnt[c] += club::color_cnt[l][c];
        vis_colors.push_back(c);
      }
    }

    for (ui a = 1; a <= k - cur_mis; ++a) {
      ui max_added = (k - cur_mis) / a;
      if (l + m + max_added <= maximum_defective.size()) {
        // pruned
        goto EXIT;
      }
      m += std::min(max_added, club::mis_cnt[a]);
      cur_mis += std::min(max_added, club::mis_cnt[a]) * a;
      assert(cur_mis <= k);

      if (cur_mis + a > k) {
        goto EXIT;
      }

      if (l == a) {
        for (ui c = 0; c < club::num_colors[0]; ++c) {
          for (ui j = 0; a + j <= k - cur_mis; ++j) {
            if (club::total_color_cnt[c] == 0) {
              m += 1;
              cur_mis += a;
              assert(cur_mis <= k);

              if (cur_mis + a > k) {
                goto EXIT;
              }
              vis_colors.push_back(c);
            } else if (a + club::total_color_cnt[c] <= k - cur_mis) {
              club::mis_cnt[a + club::total_color_cnt[c]] += 1;
            } else {
              break;
            }
            club::total_color_cnt[c] += 1;
          }
        }

        for (ui aa = a + 1; aa <= k - cur_mis; ++aa) {
          ui max_added = (k - cur_mis) / aa;
          if (l + m + max_added <= maximum_defective.size()) {
            // pruned
            goto EXIT;
          }
          m += std::min(max_added, club::mis_cnt[aa]);
          cur_mis += std::min(max_added, club::mis_cnt[aa]) * aa;
          assert(cur_mis <= k);

          if (cur_mis + aa > k) {
            goto EXIT;
          }
        }
        goto EXIT;
      } else {
        for (ui j = 0; j < h->cap(l - a); ++j) {
          ui u = h->get(l - a, j);
          ui uc = g->get_static_color(u);
          if (club::total_color_cnt[uc] == 0) {
            m += 1;
            cur_mis += a;
            assert(cur_mis <= k);

            if (cur_mis + a > k) {
              goto EXIT;
            }
            vis_colors.push_back(uc);
          } else if (a + club::total_color_cnt[uc] <= k - cur_mis) {
            club::mis_cnt[a + club::total_color_cnt[uc]] += 1;
          }
          club::total_color_cnt[uc] += 1;
        }
      }
    }
  EXIT:
    for (ui j = 2; j <= k - mis; ++j) {
      club::mis_cnt[j] = 0;
    }
    for (ui c : vis_colors) {
      club::total_color_cnt[c] = 0;
    }
    vis_colors.resize(0);
    assert(club::mis_cnt[1] == 0);

    return m;
  }
}

static ui num_maximal_defec = 0;

void move_up(ui u, const ui* updated, const ui sz) {
  if (g->neg_idx[u] > 0) {
    h->inc_all();
    for (ui i = 0; i < sz; ++i) {
      ui v = updated[i];
      h->dec(v);
    }
  } else {
    for (ui i = 0; i < sz; ++i) {
      ui v = updated[i];
      h->inc(v);
    }
  }
}

void move_down(ui u, const ui* updated, const ui sz) {
  if (g->neg_idx[u] > 0) {
    for (ui i = 0; i < sz; ++i) {
      ui v = updated[i];
      h->inc(v);
    }
    h->dec_all();
  } else {
    for (ui i = 0; i < sz; ++i) {
      ui v = updated[i];
      h->dec(v);
    }
  }
}

void reduce(ui* reduced, ui& sz) {
  const ui l = C.size() + 1;
  if (maximum_defective.size() + 1 >= k + 2 &&
      l + h->cap(l) <= maximum_defective.size()) {
    sz = 0;
    bool pruned = false;
    for (ui i = l; i-- > 0 && l - i <= k - mis;) {
      if (h->cap(i) > 0) {
        mis += l - i;
        assert(mis <= k);
        if (pruned || l + 1 + max_mis(l) <= maximum_defective.size()) {
          pruned = true;
          while (h->cap(i) > 0) {
            ui v = h->get(i, 0);
            h->pop(v);
            reduced[sz++] = v;
          }
        }
        mis -= l - i;
      }
    }
  }
}

void reduce_restore(const ui* reduced, const ui sz) {
  for (ui i = sz; i-- > 0;) {
    ui v = reduced[i];
    h->push(v, h->lab(v));
  }
}

void compute_branch(ui num_cand, const vector<ui>& nX0, vector<ui>& nbranch,
                    ui pivot) {
  const ui l = C.size();
  nbranch.resize(0);

#ifndef LISTING_MAXIMAL
#ifdef REDUCING
  if (pivot != INVALID && num_cand - prank[pivot] > h->cap(l)) {
    pivot = INVALID;
  }
#else
  pivot = INVALID;
#endif  // REDUCING
#endif  // LISTING_MAXIMAL

  if (pivot != INVALID) {
    // pivoting
    ui bidx = num_cand - prank[pivot] - 1, fidx = 0;
    nbranch.resize(num_cand - prank[pivot]);
    if (g->neg_idx[pivot] > 0) {
      for (ui j = g->cd[degidx(pivot)][pivot];
           j < g->cd[degidx(pivot)][pivot] + g->deg[degidx(pivot)][pivot];
           ++j) {
        const ui u = g->adj[j];
        if (h->lab(u) != l) {
          nbranch[bidx--] = u;
        } else {
          nbranch[fidx++] = u;
        }
      }
    } else {
      for (ui i = l + 1; i-- > 0 && l - i <= k - mis;) {
        for (ui j = 0; j < h->cap(i); ++j) {
          ui u = h->get(i, j);
          if (u == pivot || !Cuckoo::isNbr(pivot, u)) {
            if (h->lab(u) != l) {
              nbranch[bidx--] = u;
            } else {
              nbranch[fidx++] = u;
            }
          }
        }
      }
    }
    std::sort(nbranch.begin(), nbranch.begin() + fidx,
              [&](ui u, ui v) { return prank[u] < prank[v]; });
    assert(nbranch.size() == num_cand - prank[pivot]);
    assert(bidx == fidx - 1);

#ifndef LISTING_MAXIMAL
    std::sort(nbranch.begin() + fidx, nbranch.end(), [&](ui u, ui v) {
      if (h->lab(u) != h->lab(v)) {
        return h->lab(u) > h->lab(v);
      } else {
        return u > v;
      }
    });

    assert(fidx > 0);
    assert(nbranch.size() >= fidx);
    if (fidx > 1) {
      if (nbranch[fidx - 1] == pivot) {
        std::swap(nbranch[fidx - 1], nbranch[fidx - 2]);
      }
      for (ui i = fidx; i < nbranch.size(); ++i) {
        assert(nbranch[i - 1] != pivot);
        nbranch[i - 1] = nbranch[i];
      }
      nbranch.pop_back();
    } else if (nbranch.size() > 1) {
      assert(nbranch[0] == pivot);
      for (ui i = 2; i < nbranch.size(); ++i) {
        nbranch[i - 1] = nbranch[i];
      }
      nbranch.pop_back();
    }
#endif  // LISTING_MAXIMAL
  } else {
    for (ui j = 0; j < h->cap(l); ++j) {
      ui v = h->get(l, j);
      assert(!h->is_popped(v));
      nbranch.push_back(v);
    }
    std::sort(nbranch.begin(), nbranch.end(), [](ui u, ui v) { return u > v; });
    ui start = h->cap(l);
    for (ui i = l; i-- > 0 && l - i <= k - mis;) {
      bool pruned = l + max_mis_defec(l, i) <= maximum_defective.size();
      if (!pruned) {
        for (ui j = 0; j < h->cap(i); ++j) {
          ui v = h->get(i, j);
          nbranch.push_back(v);
        }
        std::sort(nbranch.begin() + start, nbranch.end(),
                  [&](ui u, ui v) { return u > v; });
        start += h->cap(i);
      } else {
        break;
      }
    }
  }
}

static bool is_skipped;

void skip_leaf(ui xlab_max, vector<ui>& skipped) {
#ifdef LISTING_MAXIMAL
#ifdef SKIP_LEAF
  const ui l = C.size();
  skipped.resize(0);
  is_skipped = false;
  if (h->cap(l) == 0 && xlab_max < l) {
    ui hlab_max = l;
    for (ui i = l; i-- > 0 && l - i <= k - mis;) {
      if (h->cap(i) > 0) {
        hlab_max = i;
        break;
      }
    }

    ui start = l > (k - mis) ? l - (k - mis) : 0;
    for (ui i = start; i < l; ++i) {
      if (mis + (l - i) + (l - xlab_max) > k &&
          mis + (l - i) + (l - hlab_max) > k) {
        while (h->cap(i) > 0) {
          ui u = h->get(i, 0);
          uint64_t prev_num_maximal = num_maximal;
          expand_clique(u, true);
          shrink_clique();
          h->pop(u);
          skipped.push_back(u);
          num_skip_leaf += 1;
          is_skipped = true;
          assert(num_maximal == prev_num_maximal + 1);
        }
      } else {
        break;
      }
    }
  }
#endif  // SKIP_LEAF
#endif  // LISTING_MAXIMAL
}

void skip_leaf_restore(const vector<ui>& skipped) {
#ifdef LISTING_MAXIMAL
#ifdef SKIP_LEAF
  const ui l = C.size();
  for (ui u : skipped) {
    h->push(u, h->lab(u));
  }
#endif  // SKIP_LEAF
#endif  // LISTING_MAXIMAL
}

static vector<vector<pair<ui, ui>>> nX;
static vector<vector<ui>> nX0;
static vector<vector<ui>> nbranch;
static vector<ui*> updated;
static vector<ui> updated_sz;
static vector<vector<ui>> skipped;

void binary_bnb(const vector<ui>& branch, vector<pair<ui, ui>>& X) {
  num_recursive_calls += 1;

  const ui l = C.size();
  bool pruned;

  ui num_proc = 0;
  ui num_cand = h->cap(l, k - mis);
  ui size_X = X.size();
  ui updated_sz;

  uint64_t prev_num_maximal = num_maximal;
  if (is_skipped) {
    prev_num_maximal -= 1;
  }

  for (ui i = 0; i < branch.size(); ++i) {
    ui u = branch[i];
    ui ul = h->lab(u);
    ui num_nbr_cand = INVALID;
    assert(mis + l - h->lab(u) <= k);

    num_cand -= 1;
    num_proc += 1;
    h->pop(u);

    if (l <= maximum_defective.size() && ul == l) {
      ui uc;
#ifdef RECOLOR
      uc = g->color[l][u];
      assert(g->color_cnt[l][uc] > 0);
      g->color_cnt[l][uc] -= 1;
      if (g->color_cnt[l][uc] == 0) {
        g->num_distinct_colors[l] -= 1;
      }
#endif  // RECOLOR
      uc = g->get_static_color(u);
      assert(club::color_cnt[l][uc] > 0);
      club::color_cnt[l][uc] -= 1;
      if (club::color_cnt[l][uc] == 0) {
        club::num_colors[l] -= 1;
      }
    }

    ui ncand1 = h->cap(l, k - mis - (l - ul));
    intersect(u, ncand1, C.size(), l - ul, updated[l], updated_sz);
    move_up(u, updated[l], updated_sz);
    if (g->neg_idx[u] > 0) {
      num_nbr_cand = ncand1 - updated_sz;
    } else {
      num_nbr_cand = updated_sz;
    }
    mis += l - ul;

    if (l + 1 > maximum_defective.size()) {
      pruned = false;
    } else {
      pruned =
          l + 1 + max_mis(l + 1) <= maximum_defective.size() ||
          l + 1 + max_mis_color(l + 1, true) <= maximum_defective.size() ||
          l + 1 + max_mis_color_dag(l + 1, true) <= maximum_defective.size();
    }
    if (!pruned) {
      ui sz2 = 0;
      ui xlab_max;
      reduce(updated[l] + updated_sz, sz2);
#ifdef LISTING_MAXIMAL
      refine(u, X, nX[l], nX0[l], xlab_max);
#endif  // LISTING_MAXIMAL

      ui ncand2 = h->cap(l + 1, k - mis);
      ui pivot;
      ui state = induce(ncand2, nX0[l], pivot);

      if (state != 1) {
        ui ncand_X = nX[l].size();
        uint64_t prev_num_maximal = num_maximal;
        expand_clique(u, ncand_X == 0 && ncand2 == 0);
        if (ncand2 > 0) {
          skip_leaf(xlab_max, skipped[l]);
          compute_branch(ncand2, nX0[l], nbranch[l], pivot);
          binary_bnb(nbranch[l], nX[l]);
          skip_leaf_restore(skipped[l]);
        } else {
          if (prev_num_maximal == num_maximal) {
            add_type1(l + 1);
          } else {
            add_type3(l + 1);
          }
          num_leaf += 1;
          num_recursive_calls += 1;
        }
        shrink_clique();
      }

      uninduce(ncand2, nX0[l]);

#ifdef LISTING_MAXIMAL
      if (num_maximal > prev_num_maximal) {
        X.emplace_back(u, ul);
      }
#endif  // LISTING_MAXIMAL

      reduce_restore(updated[l] + updated_sz, sz2);
    } else {
      add_pruned(l + 1);
    }
    mis -= l - ul;
    move_down(u, updated[l], updated_sz);
    if (g->neg_idx[u] > 0) {
      g->neg_idx[u] -= 1;
    }

    if (ul == l && num_cand == num_nbr_cand) {
      break;
    }

    if (l > maximum_defective.size()) {
      pruned = false;
    } else {
      pruned = l + max_mis(l) <= maximum_defective.size() ||
               l + max_mis_color(l, false) <= maximum_defective.size() ||
               l + max_mis_color_dag(l, false) <= maximum_defective.size();
    }

    if (pruned) {
      break;
    }
  }

  for (ui i = 0; i < num_proc; ++i) {
    ui u = branch[i];
    h->push(u, h->lab(u));
  }

  X.resize(size_X);

  if (prev_num_maximal == num_maximal) {
    add_type1(l);
  } else {
    add_type2(l);
  }
}

void initRoot(ui v, vector<pair<ui, ui>>& X) {
  static vector<bool> vis(g->n);
#ifndef NO_COLOR_REDUCTION
  static vector<bool> vis_color(g->n);
#endif  // NO_COLOR_REDUCTION

  vis[v] = true;

  vector<vector<ui>> temp(2);
  ui p = 0, q = 1;
  temp[p].reserve(g->outdeg[0][v]);
  temp[q].reserve(g->outdeg[0][v]);
  for (ui i = g->outcd[v]; i < g->outcd[v] + g->outdeg[0][v]; ++i) {
    ui u = g->outadj[i];
    vis[u] = true;
    temp[p].emplace_back(u);
  }

  bool it = true;
  while (it) {
    vector<ui>& atemp = temp[p];
    vector<ui>& btemp = temp[q];
    btemp.clear();
    it = false;
    for (auto u : atemp) {
      int d = 0;
      for (auto w : atemp) {
        if (u != w && Cuckoo::isNbr(u, w)) {
#ifndef NO_COLOR_REDUCTION
          if (!vis_color[g->color[1][w]]) {
            d++;
            vis_color[g->color[1][w]] = true;
          }
#else
          d++;
#endif  // NO_COLOR_REDUCTION
        }
      }
      if (d + k + 2 > maximum_defective.size()) btemp.emplace_back(u);
#ifndef NO_COLOR_REDUCTION
      for (auto w : atemp) {
        vis_color[g->color[1][w]] = false;
      }
#endif  // NO_COLOR_REDUCTION
    }
    p = q;
    q = 1 - p;
    if (btemp.size() != atemp.size()) it = true;
  }

  for (auto u : temp[p]) {
    h->push(u, 1);
  }

  for (ui i = g->cd[degidx(v)][v];
       i < g->cd[degidx(v)][v] + g->deg[degidx(v)][v]; ++i) {
    int u = g->adj[i];
#ifdef LISTING_MAXIMAL
    if (!vis[u]) {
      int ud = 0;
      for (auto w : temp[p]) {
        if (Cuckoo::isNbr(u, w)) {
#ifndef NO_COLOR_REDUCTION
          if (!vis_color[g->color[1][w]]) {
            ud++;
            vis_color[g->color[1][w]] = true;
          }
#else
          ud++;
#endif  // NO_COLOR_REDUCTION
        }
      }
      if (ud + k + 2 > maximum_defective.size()) {
        X.emplace_back(u, 1);
      }
#ifndef NO_COLOR_REDUCTION
      for (auto w : temp[p]) {
        vis_color[g->color[1][w]] = false;
      }
#endif  // NO_COLOR_REDUCTION
    }
#endif  // LISTING_MAXIMAL
    vis[u] = true;
  }
  temp[q].clear();

  for (ui i = 0; i < h->cap(1); ++i) {
    ui u = h->get(1, i);
    for (ui j = g->cd[degidx(u)][u];
         j < g->cd[degidx(u)][u] + g->deg[degidx(u)][u]; ++j) {
      ui w = g->adj[j];
      if (!vis[w]) {
        vis[w] = true;
        temp[q].emplace_back(w);
        int xd = 0;
        for (auto x : temp[p]) {
          if (Cuckoo::isNbr(x, w)) {
#ifndef NO_COLOR_REDUCTION
            if (!vis_color[g->color[1][x]]) {
              xd++;
              vis_color[g->color[1][x]] = true;
            }
#else
            xd++;
#endif  // NO_COLOR_REDUCTION
          }
        }
        if (xd + k + 1 > maximum_defective.size()) {
          if (w > v) {
#ifdef LISTING_MAXIMAL
            X.emplace_back(w, 0);
#endif  // LISTING_MAXIMAL
          } else {
            h->push(w, 0);
          }
        }
#ifndef NO_COLOR_REDUCTION
        for (auto x : temp[p]) {
          vis_color[g->color[1][x]] = false;
        }
#endif  // NO_COLOR_REDUCTION
      }
    }
  }
  for (auto u : temp[q]) {
    vis[u] = false;
  }
  for (ui i = g->cd[degidx(v)][v];
       i < g->cd[degidx(v)][v] + g->deg[degidx(v)][v]; ++i) {
    ui u = g->adj[i];
    vis[u] = false;
  }
  vis[v] = false;
}

void binary_bnb_outer() {
  vector<pair<ui, ui>> X;
  vector<ui> X0;
  vector<ui> nbranch;
  bool pruned;
#ifdef LISTING_MAXIMAL
  X0.reserve(g->n);
#endif  // LISTING_MAXIMAL
  nbranch.reserve(g->n);
  for (ui u = g->n; u-- > 0;) {
    if (u < maximum_defective.size()) {
      break;
    }
#ifndef NO_COLOR_REDUCING
    for (ui i = g->outcd[u]; i < g->outcd[u] + g->outdeg[0][u]; ++i) {
      ui u = g->outadj[i];
      h->push(u, 1);
    }
    induce_and_coloring_dag(1);
    while (h->cap(1) > 0) {
      h->pop(h->get(1, 0));
    }
#endif  // NO_COLOR_REDUCING
    initRoot(u, X);

    pruned = 1 + max_mis(1) <= maximum_defective.size() ||
             1 + max_mis_color(1, true) <= maximum_defective.size() ||
             1 + max_mis_color_dag(1, true) <= maximum_defective.size();

    if (!pruned) {
#ifdef LISTING_MAXIMAL
      X0.resize(0);
      for (auto [v, vl] : X) {
        if (vl == 1) {
          X0.push_back(v);
        }
      }
#endif  // LISTING_MAXIMAL

      ui ncand = h->cap(1, k);
      ui pivot;
      ui state = induce(ncand, X0, pivot);

      if (state != 1) {
        expand_clique(u, false);
        compute_branch(ncand, X0, nbranch, pivot);
        binary_bnb(nbranch, X);
        shrink_clique();
      }

      uninduce(ncand, X0);
    } else {
      add_pruned(1);
    }
    for (ui i = 0; i <= 1; ++i) {
      while (h->cap(i) > 0) {
        ui w = h->get(i, 0);
        h->pop(w);
      }
    }
    X.resize(0);
  }
}

void binary_bnb_outer_clique() {
  vector<pair<ui, ui>> X;
  vector<ui> X0;
  vector<ui> nbranch;
  ui updated_sz;
  bool pruned;

  nbranch.reserve(g->n);
  for (ui u = 0; u < g->n; ++u) {
    h->push(u, 0);
  }

  for (ui u = g->n; u-- > 0;) {
    ui ul = h->lab(u);
    ui num_nbr_cand = INVALID;

    h->pop(u);

    if (u < maximum_defective.size()) {
      continue;
    }

    ui ncand1 = h->cap(0, k);
    ui sz2;
    intersect(u, ncand1, C.size(), 0, updated[0], updated_sz);
    move_up(u, updated[0], updated_sz);

    if (g->neg_idx[u] > 0) {
      num_nbr_cand = ncand1 - updated_sz;
    } else {
      num_nbr_cand = updated_sz;
    }

    pruned = 1 + max_mis(1) <= maximum_defective.size() ||
             1 + max_mis_color(1, true) <= maximum_defective.size() ||
             1 + max_mis_color_dag(1, true) <= maximum_defective.size();

    if (!pruned) {
      ui ncand2 = h->cap(1, k - mis);
      ui pivot;
      ui state = induce(ncand2, X0, pivot);

      if (state != 1) {
        ui ncand_X = X.size();
        expand_clique(u, ncand_X == 0 && ncand2 == 0);
        compute_branch(ncand2, X0, nbranch, pivot);
        binary_bnb(nbranch, X);
        shrink_clique();
      }

      uninduce(ncand2, X0);
    } else {
      add_pruned(2);
    }
    move_down(u, updated[0], updated_sz);
    if (g->neg_idx[u] > 0) {
      g->neg_idx[u] -= 1;
    }
  }
}

vector<ui> kclique_defective(DAG* g_, ui* max_defective, ui& defective_sz,
                             ui k_, ui upper_bound) {
  std::cout << std::endl
            << "---------- Branch and bound ----------" << std::endl;
  g = g_;
  k = k_;
  max_depth = g->num_colors_static + k + 2;
  h = new BucketHeap(g->n, max_depth);
  Cuckoo::initHash(g);
  prank.resize(g->n);
  deg0.resize(g->n);

  num_removed.reserve(max_depth);

  club::color_cnt.resize(max_depth);
  club::mis_cnt.resize(k + 1);
  club::total_color_cnt.resize(max_depth);
  club::num_colors.resize(max_depth);
  club::colors.resize(max_depth);

  nX.resize(max_depth);
  nX0.resize(max_depth);
  nbranch.resize(max_depth);
  updated.resize(max_depth);
  skipped.resize(max_depth);

  srand(2);

#ifdef PRINT_STAT
  vec_type1.resize(5000);
  vec_type2.resize(5000);
  vec_type3.resize(5000);
  vec_solution.resize(5000);
  vec_pruned.resize(5000);
#endif  // PRINT_STAT

  for (ui i = 0; i < max_depth; ++i) {
    club::color_cnt[i].resize(max_depth);
#ifdef LISTING_MAXIMAL
    nX[i].reserve(g->n);
    nX0[i].reserve(g->n);
#endif  // LISTING_MAXIMAL
    nbranch[i].reserve(g->n);
    skipped[i].reserve(g->n);
    updated[i] = new ui[g->n * 2];
  }
  club::num_colors[0] = g->num_colors_static;

#ifndef LISTING_MAXIMAL
#ifndef REDUCING
  inducing = false;
#endif  // REDUCING
  for (ui i = 0; i < defective_sz; ++i) {
    maximum_defective.push_back(max_defective[i]);
  }
  // binary_bnb_outer_clique();
  if (maximum_defective.size() > k) {
    binary_bnb_outer();
  } else {
    binary_bnb_outer_clique();
  }
#else
  maximum_defective.resize(defective_sz);
#ifdef TWOHOP_REDUCTION
  binary_bnb_outer();
#else
  binary_bnb_outer_clique();
#endif  // TWOHOP_REDUCTION
#endif  // LISTING_MAXIMAL

  std::cout << "# recursive calls: "
            << Utility::integer_to_string(num_recursive_calls + num_pruned +
                                          num_skip_leaf)
            << std::endl;

#ifdef PRINT_STAT
  ui max_dep = 0;
  for (ui i = 1; i < 5000; ++i) {
    if (vec_type1[i] || vec_type2[i] || vec_type3[i] || vec_pruned[i]) {
      max_dep = i + 1;
    }
  }
  printf("# type 1 per depth: ");
  for (ui i = 1; i < max_dep; ++i) {
    printf("%lld ", vec_type1[i]);
  }
  printf("\n");
  printf("# type 2 per depth: ");
  for (ui i = 1; i < max_dep; ++i) {
    printf("%lld ", vec_type2[i]);
  }
  printf("\n");
  printf("# type 3 per depth: ");
  for (ui i = 1; i < max_dep; ++i) {
    printf("%lld ", vec_type3[i]);
  }
  printf("\n");
  printf("# solution per depth: ");
  for (ui i = 1; i < max_dep; ++i) {
    printf("%lld ", vec_solution[i]);
  }
  printf("\n");
  printf("# pruned per depth: ");
  for (ui i = 1; i < max_dep; ++i) {
    printf("%lld ", vec_pruned[i]);
  }
  printf("\n");
#endif  // PRINT_STAT

#ifdef LISTING_MAXIMAL
  std::cout << "# maximal cliques of size >= "
            << Utility::integer_to_string(defective_sz + 1) << ": "
            << Utility::integer_to_string(num_maximal) << std::endl;
#else
  std::cout << "Size of the maximum " << k
            << "-defective clique: " << maximum_defective.size() << std::endl;
#endif  // LISTING_MAXIMAL

  delete h;
  for (ui i = 0; i < max_depth; ++i) {
    delete[] updated[i];
  }
  return maximum_defective;
}

#endif  // LISTING_H_

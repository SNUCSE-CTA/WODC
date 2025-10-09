#ifndef PREPROCESS_H_
#define PREPROCESS_H_

#include "graph.h"
#include "linear_heap.h"

namespace gops {
inline void coreDecompose(const Graph *g, ui *core, ui *order);
inline void computeSuffixDefective(const Graph *g, ui k, const ui *order,
                                   ui *max_defective, ui &sz);
inline void greedyColoring(const Graph *g, const ui *order, ui *color,
                           ui &num_colors);
inline void greedyColoring(const Graph *g, ui *color, ui &num_colors);
inline void greedyColoring1(const Graph *g, const ui *order, ui *color,
                            ui &num_colors);
inline void greedyColoring1(const Graph *g, ui *color, ui &num_colors);
inline void greedyColoring2(const Graph *g, ui *color, ui &num_colors);
inline void greedyColoring2(const Graph *g, const ui *order, ui *color,
                            ui &num_colors);
inline void colorCoreDecompose(const Graph *g, const ui *color,
                               const ui num_colors, ui *core, ui *order);
}  // namespace gops

inline void gops::greedyColoring(const Graph *g, const ui *order, ui *color,
                                 ui &num_colors) {
#ifdef TIGHT_COLORING
  gops::greedyColoring2(g, order, color, num_colors);
#else
  gops::greedyColoring1(g, order, color, num_colors);
#endif  // TIGHT_COLORING
}
inline void gops::greedyColoring(const Graph *g, ui *color, ui &num_colors) {
#ifdef TIGHT_COLORING
  gops::greedyColoring2(g, color, num_colors);
#else
  gops::greedyColoring1(g, color, num_colors);
#endif  // TIGHT_COLORING
}

inline void gops::coreDecompose(const Graph *g, ui *core, ui *order) {
  ui n = g->getNumVertices();

  ui *vis = new ui[n];
  ui *degree = new ui[n];

  std::fill(vis, vis + n, 0);
  for (ui u = 0; u < n; ++u) {
    degree[u] = g->getDegree(u);
    order[u] = u;
  }

  ui max_core = 0;

  if (n > 0) {
    ListLinearHeap *heap = new ListLinearHeap(n, n - 1);
    heap->init(n, n - 1, order, degree);
    for (ui i = 0; i < n; i++) {
      ui u = 0, key = 0;
      heap->pop_min(u, key);
      if (key > max_core) max_core = key;
      core[u] = max_core;
      order[i] = u;
      if (key + i + 1 == n) {
        ui x_size = i + 1;
        heap->get_ids(order, x_size);
        assert(x_size == n);
        for (ui j = i; j < n; j++) {
          core[order[j]] = max_core;
        }
        break;
      }
      vis[u] = 1;

      for (ui j = 0; j < g->getDegree(u); j++) {
        ui v = g->getNbr(u, j);
        if (vis[v] == 0) {
          heap->decrement(v, 1);
        }
      }
    }
    delete heap;
  }

  delete[] vis;
  delete[] degree;
}

inline void gops::greedyColoring1(const Graph *g, const ui *order, ui *color,
                                  ui &num_colors) {
  const ui n = g->getNumVertices();
  std::fill(color, color + n, INVALID);

  bool *color_vis = new bool[num_colors + 1];
  std::fill(color_vis, color_vis + num_colors + 1, false);

  num_colors = 0;

  for (ui i = n; i-- > 0;) {
    ui u = order[i];
    for (ui j = 0; j < g->getDegree(u); ++j) {
      ui v = g->getNbr(u, j);
      if (color[v] == INVALID) continue;
      color_vis[color[v]] = true;
    }

    for (ui c = 0;; ++c) {
      if (!color_vis[c]) {
        color[u] = c;
        num_colors = std::max(num_colors, c);
        break;
      }
    }

    for (ui j = 0; j < g->getDegree(u); ++j) {
      ui v = g->getNbr(u, j);
      if (color[v] == INVALID) continue;
      color_vis[color[v]] = false;
    }
  }

  num_colors += 1;

  delete[] color_vis;
}

inline void gops::greedyColoring1(const Graph *g, ui *color, ui &num_colors) {
  const ui n = g->getNumVertices();
  std::fill(color, color + n, INVALID);

  bool *color_vis = new bool[num_colors + 1];
  std::fill(color_vis, color_vis + num_colors + 1, false);

  num_colors = 0;

  for (ui u = n; u-- > 0;) {
    for (ui j = 0; j < g->getDegree(u); ++j) {
      ui v = g->getNbr(u, j);
      if (color[v] == INVALID) continue;
      color_vis[color[v]] = true;
    }

    for (ui c = 0;; ++c) {
      if (!color_vis[c]) {
        color[u] = c;
        num_colors = std::max(num_colors, c);
        break;
      }
    }

    for (ui j = 0; j < g->getDegree(u); ++j) {
      ui v = g->getNbr(u, j);
      if (color[v] == INVALID) continue;
      color_vis[color[v]] = false;
    }
  }

  num_colors += 1;

  delete[] color_vis;
}

inline void gops::greedyColoring2(const Graph *g, ui *color, ui &num_colors) {
  gops::greedyColoring1(g, color, num_colors);
  ui prev_num_colors;
  ui *color_vis1 = new ui[num_colors + 1];
  ui *color_vis2 = new ui[num_colors + 1];
  ui *anchor = new ui[num_colors + 1];
  std::fill(color_vis1, color_vis1 + num_colors + 1, 0);
  std::fill(color_vis2, color_vis2 + num_colors + 1, 0);
  do {
    const ui n = g->getNumVertices();
    std::fill(color, color + n, INVALID);

    ui threshold = num_colors - 2;
    prev_num_colors = num_colors;
    num_colors = 0;

    for (ui u = n; u-- > 0;) {
      for (ui j = 0; j < g->getDegree(u); ++j) {
        ui v = g->getNbr(u, j);
        if (color[v] == INVALID) continue;
        color_vis1[color[v]] += 1;
        anchor[color[v]] = v;
      }

      for (ui c = 0;; ++c) {
        if (!color_vis1[c]) {
          color[u] = c;
          break;
        }
      }

      if (color[u] > threshold) {
        bool updated = 0;
        for (ui c = 0; c < threshold; ++c) {
          if (color_vis1[c] == 1) {
            const ui v = anchor[c];
            for (ui j = 0; j < g->getDegree(v); ++j) {
              const ui w = g->getNbr(v, j);
              if (color[w] != -1) color_vis2[color[w]] = true;
            }

            for (ui cc = 0; cc < threshold; ++cc) {
              if (c != cc && color_vis2[cc] == 0) {
                color[u] = c;
                color[v] = cc;
                color_vis1[c] = color_vis1[cc] = 0;
                color_vis2[c] = color_vis2[cc] = 0;
                updated = true;
                break;
              }
            }

            for (ui j = 0; j < g->getDegree(v); ++j) {
              const ui w = g->getNbr(v, j);
              if (color[w] != -1) color_vis2[color[w]] = 0;
            }

            if (updated) {
              break;
            }
          }
        }
      }

      num_colors = std::max(num_colors, color[u]);

      for (ui j = 0; j < g->getDegree(u); ++j) {
        ui v = g->getNbr(u, j);
        if (color[v] == INVALID) continue;
        color_vis1[color[v]] = 0;
      }
    }

    num_colors += 1;
  } while (num_colors < prev_num_colors);
  delete[] color_vis1;
  delete[] color_vis2;
  delete[] anchor;
}

inline void gops::greedyColoring2(const Graph *g, const ui *order, ui *color,
                                  ui &num_colors) {
  gops::greedyColoring1(g, order, color, num_colors);
  ui prev_num_colors;
  ui *color_vis1 = new ui[num_colors + 1];
  ui *color_vis2 = new ui[num_colors + 1];
  ui *anchor = new ui[num_colors + 1];
  std::fill(color_vis1, color_vis1 + num_colors + 1, 0);
  std::fill(color_vis2, color_vis2 + num_colors + 1, 0);
  do {
    const ui n = g->getNumVertices();
    std::fill(color, color + n, INVALID);

    ui threshold = num_colors - 2;
    prev_num_colors = num_colors;
    num_colors = 0;

    for (ui i = n; i-- > 0;) {
      ui u = order[i];
      for (ui j = 0; j < g->getDegree(u); ++j) {
        ui v = g->getNbr(u, j);
        if (color[v] == INVALID) continue;
        color_vis1[color[v]] += 1;
        anchor[color[v]] = v;
      }

      for (ui c = 0;; ++c) {
        if (!color_vis1[c]) {
          color[u] = c;
          break;
        }
      }

      if (color[u] > threshold) {
        bool updated = 0;
        for (ui c = 0; c < threshold; ++c) {
          if (color_vis1[c] == 1) {
            const ui v = anchor[c];
            for (ui j = 0; j < g->getDegree(v); ++j) {
              const ui w = g->getNbr(v, j);
              if (color[w] != -1) color_vis2[color[w]] = true;
            }

            for (ui cc = 0; cc < threshold; ++cc) {
              if (c != cc && color_vis2[cc] == 0) {
                color[u] = c;
                color[v] = cc;
                color_vis1[c] = color_vis1[cc] = 0;
                color_vis2[c] = color_vis2[cc] = 0;
                updated = true;
                break;
              }
            }

            for (ui j = 0; j < g->getDegree(v); ++j) {
              const ui w = g->getNbr(v, j);
              if (color[w] != -1) color_vis2[color[w]] = 0;
            }

            if (updated) {
              break;
            }
          }
        }
      }

      num_colors = std::max(num_colors, color[u]);

      for (ui j = 0; j < g->getDegree(u); ++j) {
        ui v = g->getNbr(u, j);
        if (color[v] == INVALID) continue;
        color_vis1[color[v]] = 0;
      }
    }

    num_colors += 1;

  } while (num_colors < prev_num_colors);
  delete[] color_vis1;
  delete[] color_vis2;
  delete[] anchor;
}

inline void gops::computeSuffixDefective(const Graph *g, ui k, const ui *order,
                                         ui *max_defective, ui &sz) {
  if (g->getNumVertices() == 0) return;

  const ui n = g->getNumVertices();
  sz = 1;

  ui *vis = new ui[n]();
  ui num_edges = 0;

  vis[order[n - 1]] = true;
  max_defective[0] = order[n - 1];

  while (sz < n) {
    ui u = order[n - sz - 1];
    for (ui i = 0; i < g->getDegree(u); ++i) {
      ui v = g->getNbr(u, i);
      if (vis[v]) {
        num_edges += 1;
      }
    }

    if (num_edges + k < (sz + 1) * sz / 2) {
      break;
    }

    vis[u] = true;
    max_defective[sz] = u;
    sz += 1;
  }

  delete[] vis;
}

inline void gops::colorCoreDecompose(const Graph *g, const ui *color,
                                     const ui num_colors, ui *core, ui *order) {
  ui n = g->getNumVertices();

  ui *vis = new ui[n]();
  ui *color_deg = new ui[n];
  ui *degree = new ui[n];
  ui *num_nbr_colors = new ui[n * num_colors]();

  for (ui u = 0; u < n; ++u) {
    ui deg = 0;
    for (ui j = 0; j < g->getDegree(u); ++j) {
      ui v = g->getNbr(u, j);
      uint64_t idx = u * num_colors + color[v];
      if (num_nbr_colors[idx] == 0) deg += 1;
      num_nbr_colors[idx] += 1;
    }
    color_deg[u] = deg;
    degree[u] = g->getDegree(u);
    order[u] = u;
  }

  ui max_color_core = 0;

  ListLinearHeap *heap = new ListLinearHeap(n, n - 1);
  ListLinearHeap *deg_heap = new ListLinearHeap(n, n - 1);
  heap->init(n, n - 1, order, color_deg);
  deg_heap->init(0, n - 1, order, degree);

  ui *vs = new ui[n];
  ui vs_size = 0;

  heap->get_min_key_ids(vs, vs_size);
  for (ui i = 0; i < vs_size; ++i) {
    ui u = vs[i];
    deg_heap->insert(u, degree[u]);
  }

  assert(vs_size > 0);
  max_color_core = color_deg[vs[0]];
  for (ui i = 0; i < n; i++) {
    ui u = 0, key = 0;
    if (deg_heap->empty()) {
      vs_size = 0;
      heap->get_min_key_ids(vs, vs_size);
      for (ui i = 0; i < vs_size; ++i) {
        ui v = vs[i];
        deg_heap->insert(v, degree[v]);
      }
      assert(vs_size > 0);
      assert(max_color_core < heap->get_key(vs[0]));
      max_color_core = heap->get_key(vs[0]);
    }
    deg_heap->pop_min(u, key, 3);
    heap->erase(u);

    assert(key == degree[u]);

    core[u] = max_color_core;
    order[i] = u;
    if (key + i + 1 == n) {
      ui x_size = i + 1;
      heap->get_ids(order, x_size);
      assert(x_size == n);
      for (ui j = i; j < n; ++j) {
        core[order[j]] = max_color_core;
      }
      break;
    }
    vis[u] = 1;

    for (ui j = 0; j < g->getDegree(u); j++) {
      ui v = g->getNbr(u, j);
      if (vis[v] == 0) {
        ui idx = v * num_colors + color[u];
        if ((--num_nbr_colors[idx]) == 0) {
          heap->decrement(v, 1);
          if (heap->get_key(v) == max_color_core) {
            deg_heap->insert(v, degree[v]);
          }
        }
        if (heap->get_key(v) <= max_color_core) {
          deg_heap->decrement(v, 1);
        }
        degree[v] -= 1;
      }
    }
  }
  delete heap;
  delete deg_heap;

  delete[] vis;
  delete[] color_deg;
  delete[] degree;
  delete[] num_nbr_colors;
  delete[] vs;
}

/**
 * from the source code of "Efficient maximum k-defective clique computation
 * with improved time complexity" (Lijun Chang, SIGMOD 2024)
 */

namespace gops {
namespace kDC {
ui max_defec_sz;
ui K;

inline void orient_graph(ui n, ui m, ui *peel_sequence, ept *pstart, ept *pend,
                         ui *edges, ui *rid) {
  for (ui i = 0; i < n; i++) rid[peel_sequence[i]] = i;
  for (ui i = 0; i < n; i++) {
    ept &end = pend[i] = pstart[i];
    for (ept j = pstart[i]; j < pstart[i + 1]; j++)
      if (rid[edges[j]] > rid[i]) edges[end++] = edges[j];
  }

#ifndef NDEBUG
  long long sum = 0;
  for (ui i = 0; i < n; i++) sum += pend[i] - pstart[i];
  assert(sum * 2 == m);
#endif
}

inline void oriented_triangle_counting(ui n, ui m, ept *pstart, ept *pend,
                                       ui *edges, ui *tri_cnt, ui *adj) {
  memset(adj, 0, sizeof(ui) * n);
  long long cnt = 0;
  memset(tri_cnt, 0, sizeof(ui) * m);
  for (ui u = 0; u < n; u++) {
    for (ept j = pstart[u]; j < pend[u]; j++) adj[edges[j]] = j + 1;

    for (ept j = pstart[u]; j < pend[u]; j++) {
      ui v = edges[j];
      for (ept k = pstart[v]; k < pend[v]; k++)
        if (adj[edges[k]]) {
          ++tri_cnt[j];
          ++tri_cnt[k];
          ++tri_cnt[adj[edges[k]] - 1];
          ++cnt;
        }
    }

    for (ept j = pstart[u]; j < pend[u]; j++) adj[edges[j]] = 0;
  }
}

bool remove_and_shrink_oriented_tri(ui &n, ui &m, ui *out_mapping,
                                    ui *peel_sequence, ept *pstart, ept *pend,
                                    ui *edges, ui *tri_cnt, ui *rid,
                                    ui *degree) {
  for (ui i = 0; i < n; i++) degree[i] = pstart[i + 1] - pstart[i];
  ept removed_edges = 0;
  for (ui i = 0; i < n; i++)
    for (ui j = pstart[i]; j < pend[i]; j++)
      if (tri_cnt[j] + K + 2 <= max_defec_sz) {
        --degree[i];
        --degree[edges[j]];
        ++removed_edges;
      }

  if (removed_edges <= m / 4) {
    return false;
  }

  ui cnt = 0;
  for (ui i = 0; i < n; i++)
    if (degree[i] > 0) {
      out_mapping[cnt] = out_mapping[i];
      rid[i] = cnt++;
    }
  ui t_cnt = 0;
  for (ui i = 0; i < n; i++)
    if (degree[peel_sequence[i]] > 0)
      peel_sequence[t_cnt++] = rid[peel_sequence[i]];
  assert(t_cnt == cnt);

#ifndef NDEBUG
  for (ui i = 0; i < n; i++) {
    ui cnt = 0;
    for (ui j = pstart[i]; j < pend[i]; j++) {
      if (tri_cnt[j] + K + 2 > max_defec_sz) ++cnt;
      assert(edges[j] < n);
    }
    assert(cnt <= degree[i]);
  }
#endif

  ui pos = 0;
  cnt = 0;
  for (ui i = 0; i < n; i++)
    if (degree[i] > 0) {
      ui start = pstart[i];
      pstart[cnt] = pos;
      for (ui j = start; j < pend[i]; j++)
        if (tri_cnt[j] + K + 2 > max_defec_sz) edges[pos++] = rid[edges[j]];
      pend[cnt] = pos;
      assert(pos - pstart[cnt] <= degree[i]);
      pos += degree[i] - (pos - pstart[cnt]);
      ++cnt;
    }
  pstart[cnt] = m = pos;
  n = cnt;

  return true;
}

void reorganize_oriented_graph(ui n, ui *tri_cnt, ept *pstart, ept *pend,
                               ept *pend2, ui *edges, ui *edges_pointer,
                               ui *buf) {
  for (ui i = 0; i < n; i++) pend2[i] = pend[i];
  for (ui i = 0; i < n; i++) {
    for (ept j = pstart[i]; j < pend[i]; j++) {
      ept &k = pend2[edges[j]];
      edges[k] = i;
      tri_cnt[k] = tri_cnt[j];
      ++k;
    }
  }

#ifndef NDEBUG
  for (ui i = 0; i < n; i++) assert(pend2[i] == pstart[i + 1]);
#endif

  for (ui i = 0; i < n; i++) {
    pend2[i] = pend[i];
    pend[i] = pstart[i];
  }
  for (ui i = 0; i < n; i++) {
    for (ept j = pend2[i]; j < pstart[i + 1]; j++) {
      ept &k = pend[edges[j]];
      edges[k] = i;
      tri_cnt[k] = tri_cnt[j];
      edges_pointer[k] = j;
      edges_pointer[j] = k;
      ++k;
    }
  }

#ifndef NDEBUG
  for (ui i = 0; i < n; i++) assert(pend[i] == pend2[i]);
#endif

  ept *ids = pend2;
  for (ui i = 0; i < n; i++) {
    if (pend[i] == pstart[i] || pend[i] == pstart[i + 1]) continue;
    ept j = pstart[i], k = pend[i], pos = 0;
    while (j < pend[i] || k < pstart[i + 1]) {
      if (k >= pstart[i + 1] || (j < pend[i] && edges[j] < edges[k])) {
        ids[pos] = edges[j];
        buf[pos++] = edges_pointer[j++];
      } else {
        ids[pos] = edges[k];
        buf[pos++] = edges_pointer[k++];
      }
    }
    assert(pos + pstart[i] == pstart[i + 1]);
    for (ept j = 0; j < pos; j++) {
      ui idx = pstart[i] + j, k = buf[j];
      edges[idx] = ids[j];
      tri_cnt[idx] = tri_cnt[k];
      edges_pointer[idx] = k;
      edges_pointer[k] = idx;
    }
  }
}

void compact_neighbors(ui u, ui *tri_cnt, ui *edges_pointer, char *deleted,
                       ept *pstart, ept *pend, ui *edges) {
  ui end = pstart[u];
  for (ui i = pstart[u]; i < pend[u]; i++)
    if (!deleted[i]) {
      edges[end] = edges[i];
      tri_cnt[end] = tri_cnt[i];
      edges_pointer[end] = edges_pointer[i];
      edges_pointer[edges_pointer[end]] = end;
      deleted[end] = 0;
      ++end;
    }
  pend[u] = end;
}

char find(ui u, ui w, ept b, ept e, char *deleted, ept &idx, ui *edges) {
  if (b >= e) return 0;

  while (b + 1 < e) {
    idx = b + (e - b) / 2;
    if (edges[idx] > w)
      e = idx;
    else
      b = idx;
  }

  idx = b;
  if (edges[idx] == w && !deleted[idx]) return 1;
  return 0;
}

void truss_peeling(ui n, ui *Qe, ui *tri_cnt, ui *edges_pointer, char *deleted,
                   ui *degree, ept *pstart, ept *pend, ui *edges) {
#ifndef NDEBUG
  char *exist = deleted;
  for (ui i = 0; i < n; i++) {
    assert(pend[i] == pstart[i + 1]);
    for (ui j = pstart[i]; j < pstart[i + 1]; j++) exist[edges[j]] = 1;
    for (ui j = pstart[i] + 1; j < pstart[i + 1]; j++)
      assert(edges[j] > edges[j - 1]);
    for (ui j = pstart[i]; j < pstart[i + 1]; j++) {
      assert(edges_pointer[edges_pointer[j]] == j);
      assert(tri_cnt[j] == tri_cnt[edges_pointer[j]]);
      assert(edges[j] != i);
      ui cnt = 0, v = edges[j];
      for (ui k = pstart[v]; k < pstart[v + 1]; k++)
        if (exist[edges[k]]) ++cnt;
      assert(cnt == tri_cnt[j]);
    }
    for (ui j = pstart[i]; j < pstart[i + 1]; j++) exist[edges[j]] = 0;
  }
  memset(deleted, 0, sizeof(char) * n);
#endif
  assert(max_defec_sz >= K + 2);
  ui t_threshold = max_defec_sz - K - 1;
  ept Qe_n = 0;
  for (ui i = 0; i < n; i++)
    for (ui j = pstart[i]; j < pend[i]; j++)
      if (tri_cnt[j] < t_threshold && edges[j] > i) {
        Qe[Qe_n++] = i;
        Qe[Qe_n++] = edges[j];
      }
  for (ept j = 0; j < Qe_n; j += 2) {
    ui u = Qe[j], v = Qe[j + 1], idx;
    find(u, v, pstart[u], pend[u], deleted, idx, edges);
    assert(edges[idx] == v);

    ui tri_n = tri_cnt[idx];
    deleted[idx] = deleted[edges_pointer[idx]] = 1;
    --degree[u];
    --degree[v];
    if (pend[u] - pstart[u] > degree[u] * 2) {
      compact_neighbors(u, tri_cnt, edges_pointer, deleted, pstart, pend,
                        edges);
      assert(degree[u] == pend[u] - pstart[u]);
    }
    if (pend[v] - pstart[v] > degree[v] * 2) {
      compact_neighbors(v, tri_cnt, edges_pointer, deleted, pstart, pend,
                        edges);
      assert(degree[v] == pend[v] - pstart[v]);
    }
    if (pend[u] - pstart[u] < pend[v] - pstart[v]) std::swap(u, v);

    if (pend[u] - pstart[u] > (pend[v] - pstart[v]) * 2) {  // binary search
      for (ept k = pstart[v]; k < pend[v]; k++)
        if (!deleted[k]) {
          if (tri_n &&
              find(u, edges[k], pstart[u], pend[u], deleted, idx, edges)) {
            --tri_n;
            --tri_cnt[edges_pointer[idx]];
            if ((tri_cnt[idx]--) == t_threshold) {
              Qe[Qe_n++] = u;
              Qe[Qe_n++] = edges[idx];
            }
            --tri_cnt[edges_pointer[k]];
            if ((tri_cnt[k]--) == t_threshold) {
              Qe[Qe_n++] = v;
              Qe[Qe_n++] = edges[k];
            }
          }
        }
    } else {  // sorted_merge
      ept ii = pstart[u], jj = pstart[v];
      while (ii < pend[u] && jj < pend[v]) {
        if (edges[ii] == edges[jj]) {
          if (!deleted[ii] && !deleted[jj]) {
            --tri_n;
            --tri_cnt[edges_pointer[ii]];
            if ((tri_cnt[ii]--) == t_threshold) {
              Qe[Qe_n++] = u;
              Qe[Qe_n++] = edges[ii];
            }
            --tri_cnt[edges_pointer[jj]];
            if ((tri_cnt[jj]--) == t_threshold) {
              Qe[Qe_n++] = v;
              Qe[Qe_n++] = edges[jj];
            }
          }

          ++ii;
          ++jj;
        } else if (edges[ii] < edges[jj])
          ++ii;
        else
          ++jj;
      }
    }
    assert(tri_n == 0);
  }
#ifndef NDEBUG
  for (ui i = 0; i < n; i++)
    for (ui j = pstart[i]; j < pend[i]; j++)
      assert(deleted[j] || tri_cnt[j] >= t_threshold);
#endif
}

inline void truss_shrink_graph(Graph *g, ui *peel_sequence, ui K_,
                               ui max_defec_sz_) {
  ui *tri_cnt = new ui[g->m];
  ui *rid = new ui[g->n];
  ui *out_mapping = new ui[g->n];
  ui *degree = new ui[g->n];
  ept *pend = new ept[g->n + 1];
  K = K_;
  max_defec_sz = max_defec_sz_;

  for (ui i = 0; i < g->n; ++i) {
    degree[i] = g->getDegree(i);
  }

  orient_graph(g->n, g->m, peel_sequence, g->pstart, pend, g->edges, rid);
  oriented_triangle_counting(g->n, g->m, g->pstart, pend, g->edges, tri_cnt,
                             rid);
  while (g->n && remove_and_shrink_oriented_tri(
                     g->n, g->m, out_mapping, peel_sequence, g->pstart, pend,
                     g->edges, tri_cnt, rid, degree)) {
    oriented_triangle_counting(g->n, g->m, g->pstart, pend, g->edges, tri_cnt,
                               rid);
  }

  ept *pend_buf = new ept[g->n + 1];
  ui *edges_pointer = new ui[g->m];
  reorganize_oriented_graph(g->n, tri_cnt, g->pstart, pend, pend_buf, g->edges,
                            edges_pointer, rid);

  ui *Qe = new ui[g->m];
  for (ui i = 0; i < g->n; i++) {
    pend[i] = g->pstart[i + 1];
    degree[i] = g->pstart[i + 1] - g->pstart[i];
  }
  char *deleted = new char[g->m];
  memset(deleted, 0, sizeof(char) * g->m);
  truss_peeling(g->n, Qe, tri_cnt, edges_pointer, deleted, degree, g->pstart,
                pend, g->edges);

  ui cnt = 0;
  for (ui i = 0; i < g->n; i++)
    if (degree[i] > 0) {
      out_mapping[cnt] = out_mapping[i];
      rid[i] = cnt++;
    }
  ui t_cnt = 0;
  for (ui i = 0; i < g->n; i++)
    if (degree[peel_sequence[i]] > 0)
      peel_sequence[t_cnt++] = rid[peel_sequence[i]];
  assert(t_cnt == cnt);
  ui pos = 0;
  cnt = 0;
  for (ui i = 0; i < g->n; i++)
    if (degree[i] > 0) {
      ui start = g->pstart[i];
      g->pstart[cnt] = pos;
      for (ui j = start; j < pend[i]; j++)
        if (!deleted[j]) {
          assert(degree[g->edges[j]] > 0);
          g->edges[pos++] = rid[g->edges[j]];
        }
      ++cnt;
    }
  g->pstart[cnt] = g->m = pos;
  g->n = cnt;

  delete[] tri_cnt;
  delete[] rid;
  delete[] out_mapping;
  delete[] degree;
  delete[] pend;
  delete[] pend_buf;
  delete[] edges_pointer;
  delete[] Qe;
  delete[] deleted;
}
}  // namespace kDC
};  // namespace gops

#endif  // PREPROCESS_H_

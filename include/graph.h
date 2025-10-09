#ifndef GRAPH_H_
#define GRAPH_H_

#include <bits/stdc++.h>

#include "utility.h"

class Graph {
 public:
  inline Graph();
  inline ~Graph();

  inline void operator=(const Graph &g) {
    n = g.n;
    m = g.m;

    pstart = new ept[n + 1];
    edges = new ui[m];

    std::copy(g.pstart, g.pstart + n + 1, pstart);
    std::copy(g.edges, g.edges + m, edges);
  }

  inline void readGraphBinary(const std::string &graph_dir);
  inline void reduce(ui lb, ui k, ui *core, ui *order);

  inline ui getNumVertices() const { return n; }
  inline ept getNumEdges() const { return m; }

  inline ui getDegree(ui u) const {
    assert(u < n);
    return pstart[u + 1] - pstart[u];
  }

  inline ui getNbr(ui u, ui i) const {
    assert(i < getDegree(u));
    return edges[pstart[u] + i];
  }

  //  private:
  ui n;
  ept m;

  ept *pstart = nullptr;
  ui *edges = nullptr;
};

inline Graph::Graph() {}

inline Graph::~Graph() {
  if (pstart != nullptr) delete[] pstart;
  if (edges != nullptr) delete[] edges;
}

inline void Graph::readGraphBinary(const std::string &graph_dir) {
  FILE *f = Utility::open_file(
      (graph_dir + std::string("/b_degree.bin")).c_str(), "rb");

  ui tt;
  fread(&tt, sizeof(int), 1, f);
  if (tt != sizeof(int)) {
    printf("sizeof int is different: edge.bin(%d), machine(%d)\n", tt,
           (int)sizeof(int));
    return;
  }
  fread(&n, sizeof(int), 1, f);
  fread(&m, sizeof(int), 1, f);

  ui *degree = new ui[n];
  fread(degree, sizeof(int), n, f);

  fclose(f);

  f = Utility::open_file((graph_dir + std::string("/b_adj.bin")).c_str(), "rb");

  if (pstart == nullptr) pstart = new ept[n + 1];
  if (edges == nullptr) edges = new ui[m];

  ui max_degree = 0;

  pstart[0] = 0;
  for (ui i = 0; i < n; i++) {
    if (degree[i] > 0) {
      fread(edges + pstart[i], sizeof(int), degree[i], f);

      // remove self loops and parallel edges
      ui *buff = edges + pstart[i];
      std::sort(buff, buff + degree[i]);
      ui idx = 0;
      for (ui j = 0; j < degree[i]; j++) {
        if (buff[j] >= n) printf("vertex id %u wrong\n", buff[j]);
        if (buff[j] == i || (j > 0 && buff[j] == buff[j - 1])) continue;
        buff[idx++] = buff[j];
      }
      degree[i] = idx;
      max_degree = std::max(max_degree, degree[i]);
    }

    pstart[i + 1] = pstart[i] + degree[i];
  }

  fclose(f);

  delete[] degree;
}

inline void Graph::reduce(ui lb, ui k, ui *core, ui *order) {
  ui nn = 0;
  ept nm = 0;
  ui min_core = 1;
  if (lb > k) {
    min_core = lb - k;
  }

  ui *degree = new ui[n]();
  ui *rank = new ui[n]();
  ui *direct_degree = new ui[n]();

  for (ui i = n; i-- > 0;) {
    ui u = order[i];
    if (core[u] < min_core) {
      break;
    }

    rank[u] = i;
    nn += 1;

    degree[u] = pstart[u + 1] - pstart[u];
    direct_degree[u] = pstart[u + 1] - pstart[u];

    for (ui j = 0; j < degree[u]; ++j) {
      ui v = edges[pstart[u] + j];
      if (rank[v]) {
        nm += 2;
        degree[v] += 1;
      } else {
        std::swap(edges[pstart[u] + j], edges[pstart[u] + degree[u] - 1]);
        degree[u] -= 1;
        direct_degree[u] -= 1;
        j -= 1;
      }
    }
  }

  ept *npstart = new ept[nn + 1];
  ui *nedges = new ui[nm];

  const ui start = n - nn;

  npstart[0] = 0;
  for (ui i = 0; i < nn; ++i) {
    ui u = order[i + start];
    npstart[i + 1] = npstart[i] + degree[u];
    rank[u] -= start;
  }

  assert(nm == npstart[nn]);

  for (ui i = 0; i < nn; ++i) {
    ui u = order[i + start];
    ui nu = rank[u];
    assert(nu == i);
    assert(nu >= 0 && nu < nn);
    for (ui j = 0; j < direct_degree[u]; ++j) {
      ui v = edges[pstart[u] + j];
      ui nv = rank[v];
      assert(nv >= 0 && nv < nn);
      nedges[npstart[nu]++] = nv;
      nedges[npstart[nv]++] = nu;
    }
  }

  for (ui i = nn; i-- > 1;) {
    npstart[i] = npstart[i - 1];
  }
  npstart[0] = 0;

  delete[] pstart;
  delete[] edges;
  delete[] degree;
  delete[] direct_degree;
  delete[] rank;

  n = nn;
  m = nm;
  pstart = npstart;
  edges = nedges;
}

#endif  // GRAPH_H_

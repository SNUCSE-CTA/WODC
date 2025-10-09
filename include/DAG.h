#ifndef DAG_H_
#define DAG_H_

#include "bucket_heap.h"
#include "graph.h"

using std::pair;
using std::vector;

class DAG {
 public:
  DAG(const Graph* g, const ui* order, ui* color_static_, const ui ncolors,
      const ui k);
  DAG(vector<pair<ui, ui>>& vset, vector<vector<ui>>& adjarr, ui e, ui ncolors,
      ui k, ui* color_static_);
  ~DAG();

  void sanity_check();

  inline ui get_static_color(const ui u);

  void resize_adj();

  ui n;   // number of nodes
  ept e;  // number of edges

  ui** outdeg;  // d[l]: (out)degrees of G_l
  ui** deg;
  ui** cd;  // cumulative degree: (starts with 0) length=n+1
  ui* outcd;
  ui* adj;  // truncated list of neighbors
  ui* outadj;
  ui* rank;  // rank of the nodes according to degeneracy ordering
  ui* idx;

  ui* neg_idx;

  ui** ssub;  // sub[l]: nodes in G_l
  ui** csub;  // sub[l]: nodes in G_l

  ui* clab;
  ui* cns;

  ui** color;
  ui** color_cnt;
  ui* num_distinct_colors;
  ui* num_colors;

  ui* color_static;
  const ui num_colors_static;

  ui* vmap;

  ui adj_off;
  ui adj_size;

  bool* inv;
  const ui k;
};

DAG::DAG(const Graph* g, const ui* order, ui* color_static_, const ui ncolors,
         const ui k)
    : num_colors_static(ncolors), k(k) {
  n = g->getNumVertices();
  e = g->getNumEdges() / 2;
  adj_off = e * 2;
  adj_size = 4 * e;

  outdeg = new ui*[ncolors + k + 2];
  deg = new ui*[ncolors + k + 2];
  cd = new ui*[ncolors + k + 2];
  outcd = new ui[n + 1]();
  adj = new ui[adj_size];
  outadj = new ui[e];
  rank = new ui[n];

  ssub = new ui*[ncolors + k + 2];
  csub = new ui*[ncolors + k + 2];
  idx = new ui[n]();

  neg_idx = new ui[n]();

  clab = new ui[n]();
  cns = new ui[ncolors + k + 2]();

  inv = new bool[n]();

  color_static = new ui[n];

  color = new ui*[ncolors + k + 2];
  color_cnt = new ui*[ncolors + k + 2];
  num_distinct_colors = new ui[ncolors + k + 2]();
  num_colors = new ui[ncolors + k + 2]();

  vmap = new ui[n];

  for (ui i = 0; i < n; ++i) {
    rank[order[i]] = n - i - 1;
    vmap[n - i - 1] = order[i];
    color_static[n - i - 1] = color_static_[vmap[n - i - 1]];
  }

  for (ui i = 1; i < ncolors + k + 2; ++i) {
    outdeg[i] = new ui[n];
    deg[i] = new ui[n];
    cd[i] = new ui[n + 1];
    ssub[i] = new ui[n];
    csub[i] = new ui[n];
    color[i] = new ui[n];
    color_cnt[i] = new ui[n]();
  }
  cns[0] = n;
  outdeg[0] = new ui[n];
  deg[0] = new ui[n];
  cd[0] = new ui[n + 1]();
  ssub[0] = new ui[n];
  csub[0] = new ui[n];
  color[0] = nullptr;

  std::iota(ssub[0], ssub[0] + n, 0);
  std::iota(csub[0], csub[0] + n, 0);

  for (ui u = 0; u < n; ++u) {
    for (ui i = 0; i < g->getDegree(u); ++i) {
      ui v = g->getNbr(u, i);
      if (rank[u] > rank[v]) {
        // make edge from u to v
        cd[0][rank[u] + 1] += 1;
        cd[0][rank[v] + 1] += 1;
        outcd[rank[u] + 1] += 1;
      }
    }
  }

  cd[0][0] = 0;
  for (ui i = 0; i < n; ++i) {
    cd[0][i + 1] += cd[0][i];
    deg[0][i] = cd[0][i + 1] - cd[0][i];
    outcd[i + 1] += outcd[i];
    outdeg[0][i] = outcd[i + 1] - outcd[i];
  }
  assert(cd[0][n] == 2 * e);
  assert(outcd[n] == e);

  for (ui u = 0; u < n; ++u) {
    for (ui i = 0; i < g->getDegree(u); ++i) {
      ui v = g->getNbr(u, i);
      adj[cd[0][rank[u]]++] = rank[v];
      if (rank[u] > rank[v]) {
        outadj[outcd[rank[u]]++] = rank[v];
      }
    }
  }

  for (ui i = n; i-- > 0;) {
    cd[0][i + 1] = cd[0][i];
    outcd[i + 1] = outcd[i];
  }
  cd[0][0] = 0;
  outcd[0] = 0;

  for (ui u = 0; u < n; ++u) {
    std::sort(adj + cd[0][u], adj + cd[0][u + 1]);
    std::sort(outadj + outcd[u], outadj + outcd[u + 1]);
  }

  for (ui i = 1; i <= k; ++i) {
    for (ui u = 0; u < n; ++u) {
      deg[i][u] = deg[0][u];
    }
  }

  sanity_check();
}

DAG::DAG(vector<pair<ui, ui>>& vset, vector<vector<ui>>& adjarr, ui e,
         ui ncolors, ui k, ui* color_static_)
    : e(e), num_colors_static(ncolors), k(k) {
  n = vset.size();
  adj_off = e * 2;
  adj_size = 4 * e;

  outdeg = new ui*[ncolors + k + 2];
  deg = new ui*[ncolors + k + 2];
  cd = new ui*[ncolors + k + 2];
  outcd = new ui[n + 1]();
  adj = new ui[adj_size];
  outadj = new ui[e];

  ssub = new ui*[ncolors + k + 2];
  csub = new ui*[ncolors + k + 2];
  idx = new ui[n]();

  neg_idx = new ui[n]();

  clab = new ui[n]();
  cns = new ui[ncolors + k + 2]();

  inv = new bool[n]();

  color = new ui*[ncolors + k + 2];
  color_cnt = new ui*[ncolors + k + 2];
  color_static = new ui[n];

  num_distinct_colors = new ui[ncolors + k + 2]();
  num_colors = new ui[ncolors + k + 2]();

  ui cap = 0;
  for (auto [v, vl] : vset) {
    if (cap < v) {
      cap = v;
    }
  }
  cap += 1;
  vmap = new ui[cap];
  rank = new ui[n];

  for (ui i = 0; i < n; ++i) {
    rank[i] = vset[i].first;
    vmap[vset[i].first] = i;
    color_static[i] = color_static_[rank[i]];
  }

  for (ui i = 1; i < ncolors + k + 2; ++i) {
    outdeg[i] = new ui[n];
    deg[i] = new ui[n];
    cd[i] = new ui[n + 1];
    ssub[i] = new ui[n];
    csub[i] = new ui[n];
    color[i] = new ui[n];
    color_cnt[i] = new ui[n]();
  }
  cns[0] = n;
  outdeg[0] = new ui[n];
  deg[0] = new ui[n];
  cd[0] = new ui[n + 1]();
  ssub[0] = new ui[n];
  csub[0] = new ui[n];
  color[0] = nullptr;

  // std::iota(sub[0], sub[0] + n, 0);
  std::iota(ssub[0], ssub[0] + n, 0);
  std::iota(csub[0], csub[0] + n, 0);

  for (ui u = 0; u < n; ++u) {
    for (ui i = 0; i < adjarr[u].size(); ++i) {
      ui v = adjarr[u][i];
      if (u > v) {
        // make edge from u to v
        cd[0][u + 1] += 1;
        cd[0][v + 1] += 1;
        outcd[u + 1] += 1;
      }
    }
  }

  cd[0][0] = 0;
  for (ui i = 0; i < n; ++i) {
    cd[0][i + 1] += cd[0][i];
    deg[0][i] = cd[0][i + 1] - cd[0][i];
    outcd[i + 1] += outcd[i];
    outdeg[0][i] = outcd[i + 1] - outcd[i];
  }
  assert(cd[0][n] == 2 * e);
  assert(outcd[n] == e);

  for (ui u = 0; u < n; ++u) {
    for (ui i = 0; i < adjarr[u].size(); ++i) {
      ui v = adjarr[u][i];
      adj[cd[0][u]++] = v;
      if (u > v) {
        outadj[outcd[u]++] = v;
      }
    }
  }

  for (ui i = n; i-- > 0;) {
    cd[0][i + 1] = cd[0][i];
    outcd[i + 1] = outcd[i];
  }
  cd[0][0] = 0;
  outcd[0] = 0;

  for (ui u = 0; u < n; ++u) {
    std::sort(adj + cd[0][u], adj + cd[0][u + 1]);
    std::sort(outadj + outcd[u], outadj + outcd[u + 1]);
  }

  for (ui i = 1; i <= k; ++i) {
    for (ui u = 0; u < n; ++u) {
      deg[i][u] = deg[0][u];
    }
  }

  sanity_check();
}

void DAG::sanity_check() {
#ifndef NDEBUG
  for (ui u = 0; u < n; ++u) {
    for (ui i = cd[0][u]; i < cd[0][u] + deg[0][u]; ++i) {
      ui v = adj[i];
      assert(color_static[u] != color_static[v]);
    }
  }
#endif  // NDEBUG
}

DAG::~DAG() {
  delete[] outcd;
  delete[] adj;
  delete[] outadj;
  delete[] rank;
  delete[] vmap;
  delete[] clab;
  delete[] cns;
  delete[] idx;
  delete[] neg_idx;

  for (ui i = 0; i < num_colors_static + k + 2; ++i) {
    delete[] outdeg[i];
    delete[] deg[i];
    delete[] cd[i];
    delete[] ssub[i];
    delete[] csub[i];
    delete[] color_cnt[i];
    if (color[i]) delete[] color[i];
  }

  delete[] outdeg;
  delete[] deg;
  delete[] cd;
  delete[] ssub;
  delete[] csub;

  delete[] color;
  delete[] color_cnt;
  delete[] num_distinct_colors;
  delete[] num_colors;

  delete[] color_static;
  delete[] inv;
}

inline ui DAG::get_static_color(const ui u) { return color_static[u]; }

void DAG::resize_adj() {
  adj_size *= 2;
  ui* tmp = adj;
  ui* nadj = new ui[adj_size];
  std::copy(tmp, tmp + adj_size / 2, nadj);
  delete[] tmp;
  adj = nadj;
}

#endif  // DAG_H_

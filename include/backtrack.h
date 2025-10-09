#ifndef BACKTRACK_H_
#define BACKTRACK_H_

#include "DAG.h"
#include "graph.h"
#include "listing.h"
#include "preprocess.h"

class Backtrack {
 public:
  Backtrack(const Graph *g, ui k);
  ~Backtrack();
  void backtrack(ui *max_defective, ui upper_bound, ui &defective_sz);
  void backtrack(ui lb);

 private:
  DAG *g;
  ui k;

  ui *core;
  ui *order;
  ui *color;

  ui num_colors;
};

Backtrack::Backtrack(const Graph *og, ui k_) {
  k = k_;
  core = new ui[og->getNumVertices()];
  order = new ui[og->getNumVertices()];
  color = new ui[og->getNumVertices()];
  gops::coreDecompose(og, core, order);
  num_colors = core[order[og->getNumVertices() - 1]];
  gops::greedyColoring(og, order, color, num_colors);

  ui degen = core[order[og->getNumVertices() - 1]];

  g = new DAG(og, order, color, num_colors, k);

#ifndef NDEBUG
  for (ui u = 0; u < g->n; ++u) {
    for (ui j = g->cd[0][u]; j < g->cd[0][u + 1]; ++j) {
      const ui v = g->adj[j];
      assert(color[g->vmap[u]] != color[g->vmap[v]]);
    }
  }
#endif  // NDEBUG
}

Backtrack::~Backtrack() {
  delete[] core;
  delete[] order;
  delete[] color;
}

void Backtrack::backtrack(ui *max_defective, ui upper_bound, ui &defective_sz) {
  auto ans = kclique_defective(g, max_defective, defective_sz, k, upper_bound);
  defective_sz = ans.size();
  for (ui i = 0; i < ans.size(); ++i) {
    max_defective[i] = ans[i];
  }
}

void Backtrack::backtrack(ui lb) {
  auto ans = kclique_defective(g, nullptr, lb, k, num_colors + k + 2);
}

#endif  // BACKTRACK_H_

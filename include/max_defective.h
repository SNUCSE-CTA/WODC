#ifndef MAX_DEFECTIVE_H_
#define MAX_DEFECTIVE_H_

#include "backtrack.h"
#include "graph.h"
#include "preprocess.h"
#include "timer.h"

namespace gops {
inline void computeMaximumDefective(const Graph* g, const ui k,
                                    ui* max_defective, ui& defective_sz,
                                    ui& upper_bound);
inline void computeMaximalDefective(const Graph* g, const ui k, const ui q);

int64_t num_backtrack = 0;
int64_t num_pruned = 0;
int64_t num_reduced = 0;

}  // namespace gops

inline void gops::computeMaximumDefective(const Graph* g_, const ui k,
                                          ui* max_defective, ui& defective_sz,
                                          ui& upper_bound) {
  Graph* g = new Graph;
  *g = *g_;

  Timer total_timer;
  Timer preproc_timer;
  ui* core = new ui[g->getNumVertices()];
  ui* order = new ui[g->getNumVertices()];
  ui* max_color_defective = nullptr;
  ui* color = nullptr;
  ui num_colors = 0;
  ui color_defective_sz = 0;
  ui ub = 0;

  ui tmp = 0;

  ui num_rounds = 1;

  ui min_truss = 0;
  ui num_edges = 0;

  Timer t;
  uint64_t elapsed;

  defective_sz = 0;
  upper_bound = g->getNumVertices();

  std::cout << "---------- Input graph ----------" << std::endl;
  std::cout << "|V|: " << Utility::integer_to_string(g->getNumVertices())
            << "; |E|: " << Utility::integer_to_string(g->getNumEdges() / 2)
            << std::endl
            << std::endl;

#ifdef CTD_MULTIROUND
  num_rounds = 3;
#endif  // CTD_MULTIROUND
#ifdef NO_PRACTICAL
  num_rounds = 0;
#endif  // NO_PRACTICAL

  max_color_defective = new ui[g->getNumVertices()];

  for (ui i = 0; i < num_rounds; ++i) {
    t.restart();
    gops::coreDecompose(g, core, order);
    elapsed = t.elapsed();

    gops::computeSuffixDefective(g, k, order, max_defective, tmp);
    defective_sz = std::max(defective_sz, tmp);
    ui max_core = core[order[g->getNumVertices() - 1]];
    ub = core[order[g->getNumVertices() - 1]] + k + 1;
    upper_bound = std::min(upper_bound, ub);
    if (defective_sz == upper_bound) {
      std::cout << "---------- Initial solution and graph reduction ----------"
                << std::endl;
      std::cout << "|S*| = " << defective_sz << std::endl;
      std::cout << "Size of the reduced graph: |V|=" << 0 << "; |E|=" << 0
                << std::endl;

      std::cout << "Size of the maximum " << k
                << "-defective clique: " << defective_sz << std::endl;
      std::cout << "Time: " << (double)preproc_timer.elapsed() / 1000000
                << " sec" << std::endl;
      goto MAXIMUM_FOUND;
    }

    g->reduce(defective_sz, k, core, order);
    if (g->getNumVertices() <= defective_sz) {
      upper_bound = defective_sz;
      std::cout << "---------- Initial solution and graph reduction ----------"
                << std::endl;
      std::cout << "|S*| = " << defective_sz << std::endl;
      std::cout << "Size of the reduced graph: |V|=" << 0 << "; |E|=" << 0
                << std::endl;

      std::cout << "Size of the maximum " << k
                << "-defective clique: " << defective_sz << std::endl;
      std::cout << "Time: " << (double)preproc_timer.elapsed() / 1000000
                << " sec" << std::endl;
      goto MAXIMUM_FOUND;
    }

    if (color == nullptr) {
      color = new ui[g->getNumVertices()];
    }
    t.restart();
    num_colors = max_core;
    gops::greedyColoring(g, color, num_colors);
    elapsed = t.elapsed();

    ub = num_colors + k;
    upper_bound = std::min(upper_bound, ub);
    if (defective_sz >= upper_bound) {
      upper_bound = defective_sz;
      std::cout << "---------- Initial solution and graph reduction ----------"
                << std::endl;
      std::cout << "|S*| = " << defective_sz << std::endl;
      std::cout << "Size of the reduced graph: |V|=" << 0 << "; |E|=" << 0
                << std::endl;
      std::cout << "Size of the maximum " << k
                << "-defective clique: " << defective_sz << std::endl;
      std::cout << "Time: " << (double)preproc_timer.elapsed() / 1000000
                << " sec"

                << std::endl;
      goto MAXIMUM_FOUND;
    }

    ui prev_num_vertices = g->getNumVertices();
    t.restart();
    gops::colorCoreDecompose(g, color, num_colors, core, order);
    elapsed = t.elapsed();
    gops::computeSuffixDefective(g, k, order, max_color_defective,
                                 color_defective_sz);
    if (color_defective_sz > defective_sz) {
      for (ui i = 0; i < color_defective_sz; ++i) {
        max_defective[i] = max_color_defective[i];
      }
      defective_sz = color_defective_sz;
    }
    ub = core[order[g->getNumVertices() - 1]] + k + 1;
    upper_bound = std::min(upper_bound, ub);
    if (defective_sz >= upper_bound) {
      upper_bound = defective_sz;
      std::cout << "---------- Initial solution and graph reduction ----------"
                << std::endl;
      std::cout << "|S*| = " << defective_sz << std::endl;
      std::cout << "Size of the reduced graph: |V|=" << 0 << "; |E|=" << 0
                << std::endl;

      std::cout << "Size of the maximum " << k
                << "-defective clique: " << defective_sz << std::endl;
      std::cout << "Time: " << (double)preproc_timer.elapsed() / 1000000
                << " sec" << std::endl;
      goto MAXIMUM_FOUND;
    }

    g->reduce(defective_sz, k, core, order);
    if (g->getNumVertices() <= defective_sz) {
      upper_bound = defective_sz;
      std::cout << "---------- Initial solution and graph reduction ----------"
                << std::endl;
      std::cout << "|S*| = " << defective_sz << std::endl;
      std::cout << "Size of the reduced graph: |V|=" << 0 << "; |E|=" << 0
                << std::endl;

      std::cout << "Size of the maximum " << k
                << "-defective clique: " << defective_sz << std::endl;
      std::cout << "Time: " << (double)preproc_timer.elapsed() / 1000000
                << " sec" << std::endl;
      goto MAXIMUM_FOUND;
    }
  }

  if (defective_sz > k + 1) {
    // truss shrink
    t.restart();
    gops::coreDecompose(g, core, order);
    gops::kDC::truss_shrink_graph(g, order, k, defective_sz);
    elapsed = t.elapsed();

    if (g->getNumVertices() == 0) {
      goto MAXIMUM_FOUND;
    }
  }

  // branch and bound
  std::cout << "---------- Initial solution and graph reduction ----------"
            << std::endl;
  std::cout << "|S*| = " << defective_sz << std::endl;
  std::cout << "Size of the reduced graph: |V|="
            << Utility::integer_to_string(g->getNumVertices())
            << "; |E|=" << Utility::integer_to_string(g->getNumEdges() / 2)
            << std::endl;
  std::cout << "Time: " << (double)preproc_timer.elapsed() / 1000000 << " sec"
            << std::endl;
#ifdef ONLY_PREPROCESSING
  exit(0);
#endif  // ONLY_PREPROCESSING
  t.restart();
  {
    Backtrack backtrack(g, k);
    backtrack.backtrack(max_defective, upper_bound, defective_sz);
  }

  elapsed = t.elapsed();

  std::cout << "Time: " << (double)elapsed / 1000000 << " sec" << std::endl;

MAXIMUM_FOUND:
  delete[] core;
  delete[] order;
  if (color) delete[] color;
  if (max_color_defective) delete[] max_color_defective;
  delete g;
  std::cout << std::endl
            << "Total Time: " << (double)total_timer.elapsed() / 1000000
            << " sec" << std::endl;
}

inline void gops::computeMaximalDefective(const Graph* g_, const ui k,
                                          const ui q) {
  Graph* g = new Graph;
  *g = *g_;

  Timer total_timer;
  Timer preproc_timer;
  ui* core = new ui[g->getNumVertices()];
  ui* order = new ui[g->getNumVertices()];
  ui* color = nullptr;
  ui num_colors = 0;
  ui ub = 0;

  ui tmp = 0;

  ui num_rounds = 1;

  Timer t;
  uint64_t elapsed;

  ui upper_bound = g->getNumVertices();

  std::cout << "---------- Input graph ----------" << std::endl;
  std::cout << "|V|: " << Utility::integer_to_string(g->getNumVertices())
            << "; |E|: " << Utility::integer_to_string(g->getNumEdges() / 2)
            << std::endl
            << std::endl;

#ifdef CTD_MULTIROUND
  num_rounds = 5;
#endif  // CTD_MULTIROUND

#ifdef NO_COLOR_REDUCTION
  num_rounds = 1;
#endif  // NO_COLOR_REDUCTION

#ifdef NO_PRACTICAL
  num_rounds = 0;
#endif  // NO_PRACTICAL
  for (ui i = 0; i < num_rounds; ++i) {
    t.restart();
    gops::coreDecompose(g, core, order);
    elapsed = t.elapsed();

    ui max_core = core[order[g->getNumVertices() - 1]];

    g->reduce(q - 1, k, core, order);

    if (g->getNumVertices() == 0) {
      break;
    }

#ifndef NO_COLOR_REDUCTION
    if (color == nullptr) {
      color = new ui[g->getNumVertices()];
    }
    t.restart();
    num_colors = max_core;
    gops::greedyColoring(g, color, num_colors);
    elapsed = t.elapsed();

    t.restart();
    gops::colorCoreDecompose(g, color, num_colors, core, order);
    elapsed = t.elapsed();

    g->reduce(q - 1, k, core, order);

    if (g->getNumVertices() == 0) {
      break;
    }
#endif  // NO_COLOR_REDUCTION
  }

#ifndef NO_COLOR_REDUCTION
#ifndef NO_PRACTICAL
  // truss shrink
  t.restart();
  gops::coreDecompose(g, core, order);
  gops::kDC::truss_shrink_graph(g, order, k, q - 1);
  elapsed = t.elapsed();

  std::cout << "---------- Graph reduction ----------" << std::endl;
  std::cout << "Size of the reduced graph: |V|="
            << Utility::integer_to_string(g->getNumVertices())
            << "; |E|=" << Utility::integer_to_string(g->getNumEdges() / 2)
            << std::endl;
  std::cout << "Time: " << (double)preproc_timer.elapsed() / 1000000 << " sec"
            << std::endl;

  if (g->getNumVertices() == 0) {
    std::cout << std::endl
              << "# maximal cliques of size >= "
              << Utility::integer_to_string(q) << ": " << 0 << std::endl;
    goto MAXIMUM_FOUND;
  }
#endif  // NO_PRACTICAL
#endif  // NO_COLOR_REDUCTION

#ifdef ONLY_PREPROCESSING
  exit(0);
#endif  // ONLY_PREPROCESSING
  t.restart();
  {
    Backtrack backtrack(g, k);
    backtrack.backtrack(q - 1);
  }

  elapsed = t.elapsed();

  std::cout << "Time: " << (double)elapsed / 1000000 << " sec" << std::endl;

MAXIMUM_FOUND:
  delete[] core;
  delete[] order;
  if (color) delete[] color;
  delete g;
  std::cout << std::endl
            << "Total Time: " << (double)total_timer.elapsed() / 1000000
            << " sec" << std::endl;
}

#endif  // MAX_DEFECTIVE_H_

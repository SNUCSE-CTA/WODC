#include "graph.h"
#include "max_defective.h"
#include "preprocess.h"

int main(int argc, char** argv) {
  Graph* g = new Graph;
  ui k = std::atoi(argv[2]);
  g->readGraphBinary(argv[1]);

#ifndef LISTING_MAXIMAL
  ui* max_defective = new ui[g->getNumVertices()];
  ui defective_sz = 0;
  ui upper_bound = 0;

#ifdef NO_REDUCTION
  gops::computeMaximumDefectiveNoReduction(g, k, max_defective, defective_sz,
                                           upper_bound);
#else
  gops::computeMaximumDefective(g, k, max_defective, defective_sz, upper_bound);
#endif  // NO_REDUCTION

  delete g;
  delete[] max_defective;
#else
  ui q = std::atoi(argv[3]);

  gops::computeMaximalDefective(g, k, q);
#endif  // LISTING_MAXIMAL
}

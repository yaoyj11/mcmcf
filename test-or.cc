#include <iostream>
#include <ctime>
#include <stdlib.h>
#include <ratio>
#include <chrono>
#include "ortools/base/commandlineflags.h"
#include "ortools/base/logging.h"
#include "ortools/graph/ebert_graph.h"
#include "ortools/graph/max_flow.h"
#include "ortools/graph/min_cost_flow.h"

namespace operations_research {

// ----- Min Cost Flow -----

// Test on a 4x4 matrix. Example taken from
// http://www.ee.oulu.fi/~mpa/matreng/eem1_2-1.htm
void MinCostFlowOn4x4Matrix() {
  LOG(INFO) << "Min Cost Flow on 4x4 Matrix";
  const int kNumSources = 4;
  const int kNumTargets = 4;
  const CostValue kCost[kNumSources][kNumTargets] = {{90, 75, 75, 80},
                                                     {35, 85, 55, 65},
                                                     {125, 95, 90, 105},
                                                     {45, 110, 95, 115}};
  const CostValue kExpectedCost = 275;
  StarGraph graph(kNumSources + kNumTargets, kNumSources * kNumTargets);
  MinCostFlow min_cost_flow(&graph);
  for (NodeIndex source = 0; source < kNumSources; ++source) {
    for (NodeIndex target = 0; target < kNumTargets; ++target) {
      ArcIndex arc = graph.AddArc(source, kNumSources + target);
      min_cost_flow.SetArcUnitCost(arc, kCost[source][target]);
      min_cost_flow.SetArcCapacity(arc, 1);
    }
  }
  for (NodeIndex source = 0; source < kNumSources; ++source) {
    min_cost_flow.SetNodeSupply(source, 1);
  }
  for (NodeIndex target = 0; target < kNumTargets; ++target) {
    min_cost_flow.SetNodeSupply(kNumSources + target, -1);
  }
  CHECK(min_cost_flow.Solve());
  CHECK_EQ(MinCostFlow::OPTIMAL, min_cost_flow.status());
  CostValue total_flow_cost = min_cost_flow.GetOptimalCost();
  CHECK_EQ(kExpectedCost, total_flow_cost);
}

}  // namespace operations_research


int main(int argc, char* argv[]) {
    gflags::ParseCommandLineFlags( &argc, &argv, true);
    operations_research::MinCostFlowOn4x4Matrix();
    return 0;
}



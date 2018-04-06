//
// Created by yaoyj11 on 4/5/18.
//
#include <map>
#include <lemon/list_graph.h>
#include "demand.h"

#ifndef MCMCF_FLOW_SOLUTION_H
#define MCMCF_FLOW_SOLUTION_H

using namespace std;
using namespace lemon;

typedef map<ListDigraph::Arc, double> Flow;
typedef ListDigraph::Arc Arc;

class FlowSolution {
public:
    Flow used_bw;
    map<Demand, Flow> flows;

    FlowSolution();

    bool empty();

    void add_flow(Demand d, Flow f);

    Flow rm_flow(Demand d);

    Flow get_flow(Demand d);

    double flow_on_edge(Arc arc);

    void update(FlowSolution sol, double theta);

};

#endif //MCMCF_FLOW_SOLUTION_H

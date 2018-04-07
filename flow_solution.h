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

typedef map<int, double> Flow;
typedef ListDigraph::Arc Arc;

class FlowSolution {
public:
    vector<double> used_bw;
    map<int, Flow> flows;

    FlowSolution(int arc_num);

    bool empty();

    void add_flow(int d, Flow f, vector<bool> &change, vector<double>&bw_change);

    Flow rm_flow(int d, vector<bool> &change, vector<double>&bw_change);

    Flow get_flow(int d);

    double flow_on_edge(int arc);

    void update(FlowSolution sol, double theta);

    void set_arc_num(int num);

};

#endif //MCMCF_FLOW_SOLUTION_H

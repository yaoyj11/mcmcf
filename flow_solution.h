//
// Created by yaoyj11 on 4/5/18.
//
#include <map>
#include<set>
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

    void add_flow(int d, Flow f, vector<double> &bw_change,  set<int>&change_edges);
    void add_flow(int d, Flow f);

    Flow rm_flow(int d, vector<double> &bw_change, set<int> &change_edges);

    Flow get_flow(int d);

    double flow_on_edge(int arc);

    void update(FlowSolution &sol, double theta);

    void set_arc_num(int num);

};

#endif //MCMCF_FLOW_SOLUTION_H

//
// Created by yaoyj11 on 4/5/18.
//

#include<string>
#include<vector>
#include<map>
#include<lemon/list_graph.h>
#include<lemon/core.h>
#include <lemon/network_simplex.h>
#include "demand.h"
#include "flow_solution.h"
#ifndef MCMCF_FRACTIONALPACKING_H
#define MCMCF_FRACTIONALPACKING_H
using namespace std;
using namespace lemon;

typedef map<Arc, int> ValMap;


class FractionalPacking {
public:
    vector<Demand> demands;
    ValMap capacity_map;
    ValMap cost_map;
    ListDigraph graph;

    FractionalPacking();

    void read_mcf(string filename);

    Flow min_cost_flow(ValMap cost, ValMap cap, int src, int dst, int val);

};


#endif //MCMCF_FRACTIONALPACKING_H

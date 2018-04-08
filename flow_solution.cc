//
// Created by yaoyj11 on 4/5/18.
//

#include "flow_solution.h"
FlowSolution::FlowSolution(int arc_num):used_bw(arc_num, 0.0),flows() {
}

bool FlowSolution::empty() {return flows.empty();}

void FlowSolution::add_flow(int x, Flow flow, vector<bool> &change, vector<double>&bw_change,
                            vector<int>&change_edges) {
    /* Fcunton to add a flow for a commodity
     * @Params:
     * x: commdofity
     * flow: arc_id -> bw
     */
    flows[x] = flow;
    for(const auto& kv :flow){
        change[kv.first]=true;
        change_edges.push_back(kv.first);
        bw_change[kv.first] += kv.second;
        used_bw[kv.first] += kv.second;
    }
}

Flow FlowSolution::rm_flow(int d, vector<bool> &change, vector<double>&bw_change, vector<int> &change_edges) {
    if(flows.count(d) ==1){
        Flow flow = flows[d];
        flows.erase(d);
        for(const auto& kv: flow){
            change[kv.first]=true;
            change_edges.push_back(kv.first);
            bw_change[kv.first] -= kv.second;
            used_bw[kv.first] -= kv.second;
        }
        return flow;
    } else{
        return Flow();
    }
}

Flow FlowSolution::get_flow(int d) {
    return flows[d];
}

double FlowSolution::flow_on_edge(int arc_id) {
    return used_bw[arc_id];
}

void FlowSolution::set_arc_num(int num) {used_bw = vector<double>(num, 0);}
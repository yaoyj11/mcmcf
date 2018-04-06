//
// Created by yaoyj11 on 4/5/18.
//

#include "flow_solution.h"
FlowSolution::FlowSolution() {
    used_bw = Flow();
    flows = map<Demand, Flow>();
}

bool FlowSolution::empty() {return flows.empty();}

void FlowSolution::add_flow(Demand x, Flow flow) {
    /* Fcunton to add a flow for a commodity
     * @Params:
     * x: commdofity
     * flow: Arc -> bw
     */
    flows[x] = flow;
    for(const auto& kv :flow){
        if (used_bw.count(kv.first) == 0){
            used_bw[kv.first] = kv.second;
        }else{
            used_bw[kv.first] += kv.second;
        }
    }
}

Flow FlowSolution::rm_flow(Demand d) {
    if(flows.count(d) ==1){
        Flow flow = flows[d];
        flows.erase(d);
        for(const auto& kv: flow){
            used_bw[kv.first] -= kv.second;
        }
        return flow;
    } else{
        return Flow();
    }
}

Flow FlowSolution::get_flow(Demand d) {
    return flows[d];
}

double FlowSolution::flow_on_edge(Arc arc) {
    return used_bw[arc];
}

void FlowSolution::update(FlowSolution sol, double theta) {
    for(const auto&kv: flows){
        Flow flow = kv.second;
        Flow new_flow = sol.flows[kv.first];
        for (const auto &edgebw: flow){
            flow[edgebw.first] = (1 - theta)* edgebw.second;
        }
        for(const auto& edgebw: new_flow){
            if(flow.count(edgebw.first) ==0){
                flow[edgebw.first] = theta* edgebw.second;
            }else{
                flow[edgebw.first] += theta* edgebw.second;
            }
        }
    }
}
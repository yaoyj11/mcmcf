//
// Created by yaoyj11 on 4/5/18.
//

#include<fstream>
#include <iostream>
#include "fractional_packing.h"

FractionalPacking::FractionalPacking() {
    capacity_map = ValMap();
    cost_map = ValMap();
    demands = vector<Demand>();
}

void FractionalPacking::read_mcf(string filename) {
    try {
        ifstream is(filename);
        int n;
        is >> n;
        std::cout<<n<<" nodes"<<endl;
        for(int i = 0; i< n; i++){
            int t;
            is>>t>>t>>t;
            graph.addNode();
        }
        is>>n;
        std::cout<<n<<" edges"<<endl;
        for(int i = 0; i< n; i++){
            int index, s, t, cap, cost;
            is>>index>>s>>t>>cap>>cost;
            ListDigraph::Arc arc = graph.addArc(graph.nodeFromId(s), graph.nodeFromId(t));
            capacity_map[arc] = cap;
            cost_map[arc] = cost;
        }
        is>>n;
        std::cout<<n<<" demands"<<endl;
        for(int i = 0; i< n; i++){
            int index, s, t, d;
            is>>index>>s>>t>>d;
            demands.push_back(Demand(s,t,d));
        }

    }catch (Exception e){
        cout<<"error reading file "<<filename<<endl;
        cout<<e.what()<<endl;
        exit(-1);
    }
}

Flow FractionalPacking::min_cost_flow(ValMap cost, ValMap cap, int src, int dst, int d) {
    NetworkSimplex<ListDigraph> ns(graph);
    ListDigraph::ArcMap<int> c(graph), u(graph);
    for (const auto &kv: capacity_map) {
        u[kv.first] = kv.second;
    }
    for (const auto &kv: cost_map) {
        c[kv.first] = kv.second;
    }
    ns.resetParams().upperMap(u).costMap(c).stSupply(graph.nodeFromId(src), graph.nodeFromId(dst), d);
    ns.run();
    Flow map;
    ListDigraph::ArcMap<int> flow_map(graph);
    ns.flowMap(flow_map);
    /*
    for (int i = 0; i < graph.maxArcId(); i++) {
        cout << i << " " << flow_map[graph.arcFromId(i)] << endl;
        if(flow_map[graph.arcFromId(i)]>0){
            cout<<"======================"<<endl;
        }
    }
     */
    return map;
}

#include <iostream>
#include <ctime>
#include <stdlib.h>
#include <ratio>
#include <chrono>
#include<lemon/list_graph.h>
#include<lemon/core.h>
#include <lemon/network_simplex.h>
#include "fractional_packing.h"

using namespace std;
using namespace std::chrono;

void test1(FractionalPacking &fp){
    int times=1000;
    std::cout<<"Testing Lemon MinCostFlow without initialization..."<<endl;
    int demand_size = 10<fp.demands.size()?10:fp.demands.size();
    for (int i = 0; i < demand_size; i++) {
        Demand d = fp.demands[i];
        //init
        ListDigraph::ArcMap<int> c(fp.graph);
        for (int j = 0; j < fp.cost_map.size(); j++) {
            c[fp.graph.arcFromId(j)] = int(fp.cost_map[j]);
        }
        fp.network_simplex->costMap(c);
        ListDigraph::ArcMap<int> u(fp.graph);
        for (int j = 0; j < fp.capacity_map.size(); j++) {
            u[fp.graph.arcFromId(j)] = fp.capacity_map[j];
        }
        fp.network_simplex->upperMap(u);
        ListDigraph::NodeMap<int> demand(fp.graph);
        demand[fp.graph.nodeFromId(d.src)] = d.val;
        demand[fp.graph.nodeFromId(d.dst)] = -d.val;
        fp.network_simplex->supplyMap(demand);
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        for(int k=0;k<times;k++){
            fp.network_simplex->run();
        }
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t2 - t1;
        std::cout << "It took me " << time_span.count() / 1000/times << " milliseconds." << endl;
    }
}

void test2(FractionalPacking &fp){
    int times=1000;
    std::cout<<"Testing Lemon MinCostFlow with initialization..."<<endl;
    int demand_size = 10<fp.demands.size()?10:fp.demands.size();
    for (int i = 0; i < demand_size; i++) {
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        Demand d = fp.demands[i];
        for(int j=0;j<times;j++){
            //init
            ListDigraph::ArcMap<int> c(fp.graph);
            for (int k = 0; k < fp.cost_map.size(); k++) {
                c[fp.graph.arcFromId(i)] = int(fp.cost_map[k]);
            }
            fp.network_simplex->costMap(c);
            ListDigraph::ArcMap<int> u(fp.graph);
            for (int k = 0; k < fp.capacity_map.size(); k++) {
                u[fp.graph.arcFromId(k)] = fp.capacity_map[k];
            }
            fp.network_simplex->upperMap(u);
            ListDigraph::NodeMap<int> demand(fp.graph);
            demand[fp.graph.nodeFromId(d.src)] = d.val;
            demand[fp.graph.nodeFromId(d.dst)] = -d.val;
            fp.network_simplex->supplyMap(demand);
            fp.network_simplex->run();
        }
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t2 - t1;
        std::cout << "It took me " << time_span.count() / 1000/times << " milliseconds." << endl;
    }

}

int main(int argc, char* argv[]) {
    srand(41);
    std::cout << "Hello, World!" << std::endl;
    string filenumber;
    if(argc>1){
        filenumber = string(argv[1]);
    }
    else {
        filenumber = "743426";
    }
    FractionalPacking fp("/Users/yaoyj11/project/mcmcf/data/test-mcf"+filenumber+".net");
    fp.network_simplex = new NetworkSimplex<ListDigraph>(fp.graph);
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    fp.time_debug=true;
    test1(fp);
    test2(fp);
    return 0;
}



#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>
#include "fractional_packing.h"

using namespace std;
using namespace std::chrono;

int main() {

    std::cout << "Hello, World!" << std::endl;
    string filenumber;
    int cost;
    //cin>>filenumber>>cost;
    filenumber="796829";
    cost = 78134;
    //lptime: 72s
    //python time 696
    FractionalPacking fp("/home/yaoyj11/project/mcmcf/data/test-mcf"+filenumber+".net");
    fp.set_buget(cost);
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    fp.fractional_packing(0.05);
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double, std::micro> time_span = t2 - t1;
    std::cout << "It took me " << time_span.count() / 1000000 << " seconds." << endl;
    cout<<"Number of min_cost_flows: "<<fp.min_cost_count<<endl;
    cout<<"Cost: "<<fp.get_cost()<<endl;

/*
for (int i = 0; i < fp.demands.size(); i++) {
    Demand d = fp.demands[i];
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    int times=100;
    for (int j = 0; j <times; j++) {
        fp.min_cost_flow(d.src, d.dst, d.val, fp.cost_map, fp.capacity_map);
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double, std::micro> time_span = t2 - t1;
    std::cout << "It took me " << time_span.count() / 1000/times << " milliseconds." << endl;
}
*/
return 0;
}

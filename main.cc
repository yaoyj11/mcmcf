#include <iostream>
#include <ctime>
#include <stdlib.h>
#include <ratio>
#include <chrono>
#include "fractional_packing.h"

using namespace std;
using namespace std::chrono;

int main(int argc, char* argv[]) {
    srand(41);
    std::cout << "Hello, World!" << std::endl;
    string filenumber;
    int cost;
    //cin>>filenumber>>cost;
    //filenumber="796829";
    //cost = 78134;
    cost = 74600;
    //lptime: 72s
    //python time 696
    //mcmcf time 13s

    //filenumber="478237";
    //cost = 84121;
    //lptime: 64.95
    //python time: 739
    //filenumber="478237";
    //cost = 84121;
    if(argc>1){
        filenumber = string(argv[1]);
    }
    else {
        filenumber = "738086";
    }
    FractionalPacking fp("/home/yaoyj11/project/mcmcf/data/test-mcf"+filenumber+".net");
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    fp.time_debug=true;
    bool suc = true;

    if(argc>2){
        cost=atoi(argv[2]);
        suc=fp.fractional_packing(cost, 0.05);
    }else{
       cost =fp.min_cost(0.05);
    }
    //cost=68597.74
    cout<<"times"<<endl;
    for(int i=0; i<fp.mab_times.size();i++){
        cout<<i<<" "<<fp.mab_times[i]<<" "<<fp.mab_index[i]<<endl;
    }
    if(fp.time_debug){

        cout<<"init flow time: "<<fp.init_flow_time<<endl;
        cout<<"min cost time: "<<fp.min_cost_time<<endl;

        cout<<"new_ton_time: "<<fp.new_ton_time<<endl;

        cout<<"compute potential_time: "<<fp.compute_potential_time<<endl;
        cout<<"update potential_time: "<<fp.update_potential_time<<endl;

        cout<<"update_flow_time: "<<fp.update_flow_time<<endl;
        cout<<"iteration_time: "<<fp.iteration_time<<endl;
        cout<<"iteration_all_time: "<<fp.iteration_all_time<<endl;
        cout<<"draw index time: "<<fp.draw_index_time<<endl;

        cout<<"succcess? "<<suc<<endl;
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double, std::micro> time_span = t2 - t1;
    std::cout << "It took me " << time_span.count() / 1000000 << " seconds." << endl;
    cout<<"Number of min_cost_flows: "<<fp.min_cost_count<<endl;
    cout<<"Number of updates: "<<fp.update_count<<endl;
    cout<<"Cost: "<<cost<<endl;
    cout<<endl;

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

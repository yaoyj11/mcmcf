#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>
#include "fractional_packing.h"
using  namespace std;
using namespace std::chrono;
int main() {

    std::cout << "Hello, World!" << std::endl;
    FractionalPacking fp;
    fp.read_mcf("/home/yaoyj11/project/mcmcf/data/test-mcf244224.net");
    for(int i = 0; i<fp.demands.size(); i++) {
        Demand d = fp.demands[i];
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        for (int j = 0; j < 100; j++) {
            fp.min_cost_flow(fp.cost_map, fp.capacity_map, d.src, d.dst, d.val);
        }
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t2 - t1;
        std::cout << "It took me " << time_span.count()/100000 << " microseconds."<<endl;
    }
    return 0;
}
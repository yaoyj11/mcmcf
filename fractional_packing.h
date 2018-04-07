//
// Created by yaoyj11 on 4/5/18.
//

#include<string>
#include<vector>
#include <time.h>
#include <ctime>
#include <ratio>
#include <chrono>
#include<map>
#include<cmath>
#include <stdexcept>
#include<lemon/list_graph.h>
#include<lemon/core.h>
#include <lemon/network_simplex.h>
#include "demand.h"
#include "edge.h"
#include "flow_solution.h"

#ifndef MCMCF_FRACTIONALPACKING_H
#define MCMCF_FRACTIONALPACKING_H
using namespace std;
using namespace lemon;


class FractionalPacking {
public:

    double budget;

    int min_cost_count;

    int update_count;

    int _i;

    int _m;

    vector<double>_u;

    vector<double>_f;

    vector<double>bw_change;

    vector<bool> change;

    double _epsilon;

    double _alpha;

    double _potential;

    double _rou;

    bool time_debug;

    double init_flow_time;

    double min_cost_time;

    double new_ton_time;

    double potential_time;

    double update_flow_time;

    double iteration_time;

    vector<Demand> demands;

    vector<Edge> edges;

    vector<int> capacity_map;

    vector<double> cost_map;

    vector<double> beta;

    ListDigraph graph;

    FlowSolution solution;

    FractionalPacking();

    FractionalPacking(string file_name);

    void read_mcf(string filename);

    Flow min_cost_flow(int src, int dst, int val, vector<double> cost, vector<int> cap);

    void set_buget(double b);

    FlowSolution fractional_packing(double epsilon, bool restart= true);

    double get_cost();

private:

    void compute_init_flow();

    string current_date_time();

    double compute_potential_function();

    double update_potential_function();

    void iteration();

    Flow update_flow(Flow oldx, Flow newx, double theta);

    double compute_theta_newton_raphson(double theta0, double mintheta, double
    maxtheta);

    double update_theta(double theta, vector<double> ax, vector<double> ax_star);

};


#endif //MCMCF_FRACTIONALPACKING_H

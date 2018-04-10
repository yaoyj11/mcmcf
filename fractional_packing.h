//
// Created by yaoyj11 on 4/5/18.
//

#include<string>
#include<vector>
#include <time.h>
#include <stdlib.h>
#include <ctime>
#include <ratio>
#include <chrono>
#include <cassert>
#include<map>
#include<set>
#include<cmath>
#include <stdexcept>
#include<lemon/list_graph.h>
#include<lemon/core.h>
#include <lemon/network_simplex.h>
#include "demand.h"
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

    vector<double>_y;

    double delta_phi_x;

    vector<double>bw_change;

    set<int> change_edges;

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

    double iteration_all_time;

    double draw_index_time;

    vector<Demand> demands;

    vector<int> capacity_map;

    vector<double> inverse_capacity;

    vector<double> cost_map;

    vector<double> beta;

    vector<bool> improve;

    bool improve_flag;

    ListDigraph graph;

    NetworkSimplex<ListDigraph> *network_simplex;

    ListDigraph::ArcMap<int> *dual_cost;

    ListDigraph::ArcMap<int> *relax_cap;

    FlowSolution solution;

    FractionalPacking();

    FractionalPacking(string file_name);

    ~FractionalPacking();

    void read_mcf(string filename);

    Flow min_cost_flow(int src, int dst, int val, const vector<double> &cost, vector<int> cap = vector<int>());

    Flow min_cost_flow(int src, int dst, int d, ListDigraph::ArcMap<int> *c=0);

    Flow min_cost_flow(int src, int dst, int d, int scale, ListDigraph::ArcMap<int> *c, ListDigraph::ArcMap<int> *
            cap);

    int min_cost_flow_cost(int src, int dst, int d);

    double min_cost_flow_cost(int src, int dst, int d, int scale, ListDigraph::ArcMap<int> *c,
                           ListDigraph::ArcMap<int> * cap);

    bool fractional_packing(double b, double epsilon, bool restart= true);

    double min_cost(double epsilon);

    double get_cost();

private:

    void set_buget(double b);

    void compute_init_flow();

    string current_date_time();

    double compute_potential_function(bool recompute_u=false);

    double update_potential_function();

    void iteration();

    bool iteration_all();

    Flow update_flow(const Flow &oldx, const Flow &newx, double theta);

    double compute_theta_newton_raphson(vector<double> & ax, vector<double> &ax_star, double theta0, double mintheta, double maxtheta);

    inline double update_theta(double theta, const vector<double> &ax, const vector<double> &ax_star);

    vector<double> mab_reward;

    vector<double> mab_index;

    double mab_average;

    bool mab_flag;

    vector<int> mab_times;

    int draw_demand_index();

    void update_mab(int index, double r);

    double compute_delta_phi_x();
};


#endif //MCMCF_FRACTIONALPACKING_H

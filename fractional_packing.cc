//
// Created by yaoyj11 on 4/5/18.
//

#include<fstream>
#include <iostream>
#include "fractional_packing.h"

using namespace std::chrono;

FractionalPacking::FractionalPacking() : budget(0), _i(0), min_cost_count(0), capacity_map(), cost_map(), demands(),
                                         beta(), solution(0) {
    _m = 0;
    _epsilon = 2.0;
    _alpha = 0;
    _potential = 0;
    _rou = 0;
    update_count = 0;
    network_simplex = 0;
    dual_cost = 0;
    relax_cap = 0;
    delta_phi_x = 0;
}

FractionalPacking::FractionalPacking(string filename) : budget(0), _i(0), min_cost_count(0), capacity_map(),
                                                        cost_map(), demands(), beta(), solution(0) {
    network_simplex = 0;
    dual_cost = 0;
    relax_cap = 0;
    update_count = 0;
    _m = 0;
    _epsilon = 2.0;
    _alpha = 0;
    _potential = 0;
    _rou = 0;
    delta_phi_x = 0;
    try {
        ifstream is(filename);
        int n;
        is >> n;
        std::cout << n << " nodes" << endl;
        for (int i = 0; i < n; i++) {
            int t;
            is >> t >> t >> t;
            graph.addNode();
        }
        is >> n;
        std::cout << n << " edges" << endl;
        solution.set_arc_num(n);
        capacity_map = vector<int>(n, 0);
        inverse_capacity = vector<double>(n, 0);
        cost_map = vector<double>(n, 0);
        for (int i = 0; i < n; i++) {
            int index, s, t, cap, cost;
            is >> index >> s >> t >> cap >> cost;
            ListDigraph::Arc arc = graph.addArc(graph.nodeFromId(s), graph.nodeFromId(t));
            capacity_map[i] = cap;
            inverse_capacity[i] = 1.0 / cap;
            cost_map[i] = cost;
        }
        is >> n;
        std::cout << n << " demands" << endl;
        for (int i = 0; i < n; i++) {
            int index, s, t, d;
            is >> index >> s >> t >> d;
            demands.push_back(Demand(s, t, d));
        }

    } catch (Exception e) {
        cout << "error reading file " << filename << endl;
        cout << e.what() << endl;
        exit(-1);
    }
}

FractionalPacking::~FractionalPacking() {
    if (network_simplex != 0) {
        delete network_simplex;
        network_simplex = 0;
    }
    if (dual_cost != 0) {
        delete dual_cost;
        dual_cost = 0;
    }
    if (relax_cap != 0) {
        delete relax_cap;
        relax_cap = 0;
    }
}

void FractionalPacking::read_mcf(string filename) {
    try {
        ifstream is(filename);
        int n;
        is >> n;
        std::cout << n << " nodes" << endl;
        for (int i = 0; i < n; i++) {
            int t;
            is >> t >> t >> t;
            graph.addNode();
        }
        is >> n;
        std::cout << n << " edges" << endl;
        capacity_map = vector<int>(n, 0);
        inverse_capacity = vector<double>(n, 0);
        cost_map = vector<double>(n, 0);
        solution.set_arc_num(n);
        for (int i = 0; i < n; i++) {
            int index, s, t, cap, cost;
            is >> index >> s >> t >> cap >> cost;
            ListDigraph::Arc arc = graph.addArc(graph.nodeFromId(s), graph.nodeFromId(t));
            capacity_map[i] = cap;
            inverse_capacity[i] = 1.0 / cap;
            cost_map[i] = cost;
        }
        is >> n;
        std::cout << n << " demands" << endl;
        for (int i = 0; i < n; i++) {
            int index, s, t, d;
            is >> index >> s >> t >> d;
            demands.push_back(Demand(s, t, d));
        }
    } catch (Exception e) {
        cout << "error reading file " << filename << endl;
        cout << e.what() << endl;
        exit(-1);
    }
}

double FractionalPacking::min_cost(double epsilon) {
    _m = cost_map.size() + 1;
    _u = vector<double>(_m, 0);
    _y = vector<double>(_m, 0);
    _f = vector<double>(_m, 0);
    _fu=vector<double>(_m,0);
    _bb = vector<double>(cost_map.size(),0);
    improve= vector<bool>(demands.size(),true);
    improve_flag=false;
    bw_change = vector<double>(cost_map.size(), 0);
    _epsilon = cost_map.size() - 1;
    network_simplex = new NetworkSimplex<ListDigraph>(graph);
    dual_cost = new ListDigraph::ArcMap<int>(graph);
    relax_cap = new ListDigraph::ArcMap<int>(graph);
    mab_reward = vector<double>(demands.size(), 0);
    mab_index = vector<double>(demands.size(), 0);
    mab_times = vector<int>(demands.size(), 0);
    mab_average = 0;
    mab_flag = true;
    compute_init_flow();
    double min = get_cost();
    double max = -1;
    //first determine a upper bound
    cout<<"try "<<2*min<<endl;
    while (fractional_packing(2 * min, epsilon, false)!=0) {
        min = 2 * min;
        cout<<"try "<<2*min<<endl;
        cout << "min: " << min << "max: " << max << endl;
    }
    FlowSolution res = solution;
    max = 2 * min * (1 + epsilon);
    cout << "min: " << min << "max: " << max << endl;
    while (max >min&& get_cost()>(min*(1+epsilon))) {
        double trial = (max + min) / 2;
        cout << "try " << trial << endl;
        int ret = fractional_packing(trial, epsilon, false);
        if (ret == 0) {
            max = trial * (1 - epsilon * epsilon);
            cout << "suc" << endl;
            res = solution;
        } else if (ret == 1) {
            min = trial;
            cout << "fail" << endl;
        } else {
            min = trial * (1 + epsilon);
            cout << "no improvement" << endl;
            cout << "min: " << min << "max: " << max << endl;
        }
    }
    solution = res;

    return get_cost();
}

int FractionalPacking::fractional_packing(double b, double epsilon, bool restart) {
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    set_buget(b);
    _epsilon = cost_map.size()-1;

    if (restart) {
        _m = cost_map.size() + 1;
        _u = vector<double>(_m, 0);
        _y = vector<double>(_m, 0);
        _fu=vector<double>(_m,0);
        _bb = vector<double>(cost_map.size(),0);
        _f = vector<double>(_m, 0);
        improve= vector<bool>(demands.size(),true);
        improve_flag=false;
        bw_change = vector<double>(cost_map.size(), 0);
        _epsilon = cost_map.size() - 1;
        network_simplex = new NetworkSimplex<ListDigraph>(graph);
        dual_cost = new ListDigraph::ArcMap<int>(graph);
        relax_cap = new ListDigraph::ArcMap<int>(graph);
        mab_reward = vector<double>(demands.size(), 0);
        mab_index = vector<double>(demands.size(), 0);
        mab_times = vector<int>(demands.size(), 0);
        mab_average = 0;
        mab_flag = true;
        compute_init_flow();
    }
    theoritical_alpha = 1 / _epsilon * log(3 * _m);
    _alpha = theoritical_alpha;
    compute_potential_function(true);
    while (_epsilon - epsilon> -1e-6) {
        _adaptive_alpha = alpha_times;
        alpha_flag = true;
        cout << current_date_time() << " epsilon: " << _epsilon << endl;
        while (_potential > _bound && _rou>(1+epsilon)) {
            //while(_rou>1+_epsilon){
            if(rand()%(3*demands.size())!=0) {
                iteration();
            }else {
                int res = iteration_all();
                if(res!=0){
                    if (time_debug) {
                        high_resolution_clock::time_point t3 = high_resolution_clock::now();
                        duration<double, std::micro> time_span = t3 - t2;
                        cout<<"fractional time" <<time_span.count() / 1000<<endl;
                    }
                    return res;
                }
            }
        }
        if (abs(_epsilon - epsilon)<1e-6) {
            compute_potential_function(true);
            cout<<"succeed: _rou "<<_rou<< " cost"<<get_cost()<<endl;
            break;
        }
        _epsilon /= 2.0;
        if (_epsilon < epsilon) {
            _epsilon = epsilon;
        }
        theoritical_alpha = 1 / _epsilon * log(3 * _m);
        _alpha = theoritical_alpha;
        compute_potential_function();
    }
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        cout<<"fractional time" <<time_span.count() / 1000<<endl;
    }
    return 0;
}

Flow FractionalPacking::min_cost_flow(int src, int dst, int d, const vector<double> &cost,
                                      vector<int> cap) {
    min_cost_count++;
    //TODO: Optimize min_cost_flow
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    ListDigraph::ArcMap<int> c(graph);
    for (int i = 0; i < cost.size(); i++) {
        c[graph.arcFromId(i)] = int(cost[i]);
    }
    network_simplex->costMap(c);
    if (!cap.empty()) {
        ListDigraph::ArcMap<int> u(graph);
        for (int i = 0; i < cap.size(); i++) {
            u[graph.arcFromId(i)] = cap[i];
        }
        network_simplex->upperMap(u);
    }
    ListDigraph::NodeMap<int> demand(graph);
    demand[graph.nodeFromId(src)] = d;
    demand[graph.nodeFromId(dst)] = -d;
    network_simplex->supplyMap(demand);
    network_simplex->run();
    Flow map;
    ListDigraph::ArcMap<int> flow_map(graph);
    network_simplex->flowMap(flow_map);
    for (int i = 0; i < graph.maxArcId(); i++) {
        if (flow_map[graph.arcFromId(i)] > 0) {
            map[i] = flow_map[graph.arcFromId(i)];
        }
    }
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        min_cost_time += time_span.count() / 1000;
    }
    return map;
}

Flow FractionalPacking::min_cost_flow(int src, int dst, int d, ListDigraph::ArcMap<int> *c) {
    min_cost_count++;
    //TODO: Optimize min_cost_flow
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    if (c != 0) {
        network_simplex->costMap(*c);
    }
    ListDigraph::NodeMap<int> demand(graph);
    demand[graph.nodeFromId(src)] = d;
    demand[graph.nodeFromId(dst)] = -d;
    network_simplex->supplyMap(demand);
    network_simplex->run();
    Flow map;
    ListDigraph::ArcMap<int> flow_map(graph);
    network_simplex->flowMap(flow_map);
    for (int i = 0; i < graph.maxArcId(); i++) {
        if (flow_map[graph.arcFromId(i)] > 0) {
            map[i] = flow_map[graph.arcFromId(i)];
        }
    }
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        min_cost_time += time_span.count() / 1000;
    }
    return map;
}

Flow FractionalPacking::min_cost_flow(int src, int dst, int d, int scale, ListDigraph::ArcMap<int> *c,
        ListDigraph::ArcMap<int> *
cap) {
    min_cost_count++;
    //TODO: Optimize min_cost_flow
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    if (c != 0) {
        network_simplex->costMap(*c);
    }
    if (cap != 0) {
        network_simplex->upperMap(*cap);
    }
    ListDigraph::NodeMap<int> demand(graph);
    demand[graph.nodeFromId(src)] = d*scale;
    demand[graph.nodeFromId(dst)] = -d*scale;
    network_simplex->supplyMap(demand);
    network_simplex->run();
    Flow map;
    ListDigraph::ArcMap<int> flow_map(graph);
    network_simplex->flowMap(flow_map);
    double inv_s=1/double(scale);
    for (int i = 0; i < graph.maxArcId(); i++) {
        if (flow_map[graph.arcFromId(i)] > 0) {
            map[i] = flow_map[graph.arcFromId(i)]*inv_s;
        }
    }
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        min_cost_time += time_span.count() / 1000;
    }
    return map;
}

int FractionalPacking::min_cost_flow_cost(int src, int dst, int d) {
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    min_cost_count++;
    //TODO: Optimize min_cost_flow
    ListDigraph::NodeMap<int> demand(graph);
    demand[graph.nodeFromId(src)] = d;
    demand[graph.nodeFromId(dst)] = -d;
    network_simplex->supplyMap(demand);
    network_simplex->run();
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        min_cost_time += time_span.count() / 1000;
    }
    return network_simplex->totalCost();
}

double FractionalPacking::min_cost_flow_cost(int src, int dst, int d, int scale, ListDigraph::ArcMap<int> *c,
                                          ListDigraph::ArcMap<int> *cap) {
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    min_cost_count++;
    //TODO: Optimize min_cost_flow
    if (c != 0) {
        network_simplex->costMap(*c);
    }
    if (cap != 0) {
        network_simplex->upperMap(*cap);
    }
    ListDigraph::NodeMap<int> demand(graph);
    demand[graph.nodeFromId(src)] = d*scale;
    demand[graph.nodeFromId(dst)] = -d*scale;
    network_simplex->supplyMap(demand);
    network_simplex->run();
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        min_cost_time += time_span.count() / 1000;
    }
    return network_simplex->totalCost()/double(scale);
}

void FractionalPacking::set_buget(double b) {
    budget = b;
    beta = vector<double>(cost_map.size(), 0);
    for (int i = 0; i < cost_map.size(); i++) {
        beta[i] = cost_map[i] / budget;
    }
}


void FractionalPacking::compute_init_flow() {
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    solution.set_arc_num(cost_map.size());
    for (int i = 0; i < demands.size(); i++) {
        Demand d = demands[i];
        Flow f = min_cost_flow(d.src, d.dst, d.val, cost_map, capacity_map);
        if (f.empty()) {
            throw std::runtime_error("No feasible solution for single commodity");
        }
        solution.add_flow(i, f, bw_change, change_edges);
    }
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        init_flow_time += time_span.count() / 1000;
    }
}

string FractionalPacking::current_date_time() {
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
    return buf;
}

double FractionalPacking::compute_potential_function(bool recompute_u) {
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    _bound = exp(_alpha * _epsilon);
    if (recompute_u) {
        _u[_u.size() - 1] = 0;
        for (int i = 0; i < cost_map.size(); i++) {
            _u[i] = solution.used_bw[i] * inverse_capacity[i];
            _bb[i]=solution.used_bw[i]*beta[i];
            _u[_u.size() - 1] += _bb[i];
            bw_change[i] = 0;
        }
    }
    change_edges.clear();
    _rou = 0;
    for (int i = 0; i < _u.size(); i++) {
        _rou = _rou > _u[i] ? _rou : _u[i];
    }
    _potential = 0;
    /*
    double sum_y = 0;
    for (int i = 0; i < _u.size(); i++) {
        sum_y += _alpha * _f[i];
    }
    sum_y += delta_phi_x/9;
     */
    delta_phi_x = 0;
    for (int i = 0; i < _u.size(); i++) {
        _f[i] = exp(_alpha * (_u[i] - 1.0));
        _fu[i]=_f[i]*_u[i];
        _potential += _f[i];
        delta_phi_x += _f[i] * _u[i];
        _y[i] = _alpha * _f[i];
    }
    _y[_y.size() - 1] += delta_phi_x/9;

    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        compute_potential_time += time_span.count() / 1000;
    }
    return _potential;
}

double FractionalPacking::update_potential_function() {
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    for (const auto &i:change_edges) {
        _u[i] = solution.used_bw[i] * inverse_capacity[i];
        _f[i] = exp(_alpha * (_u[i] - 1.0));
        _fu[i]=_f[i]*_u[i];
        _y[i] = _alpha * _f[i];
        _bb[i] = solution.used_bw[i] * beta[i];
        bw_change[i] = 0;
    }
    _u[_m-1]=0;
    for(int i=0; i<cost_map.size();i++){
        _u[_m-1]+=_bb[i];
    }
    _f[_f.size() - 1] = exp(_alpha * (_u.back() - 1.0));
    _fu[_m-1]=_f[_m-1]*_u[_m-1];
    _y[_m-1] = _alpha*_f[_m-1];
    _rou = 0;

    for (int i = 0; i < _u.size(); i++) {
        _rou = _rou > _u[i] ? _rou : _u[i];
    }
    _potential = 0;
    delta_phi_x=0;
    for(int i=0; i<_m; i++){
        delta_phi_x+=_fu[i];
        _potential +=_f[i];
    }
    _y[_y.size() - 1] += delta_phi_x/9;
    change_edges.clear();
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        update_potential_time += time_span.count() / 1000;
    }
    return _potential;
}

double FractionalPacking::compute_delta_phi_x() {
    double delta_phi_x = 0;
    //compute new cost in (4)
    for (int i = 0; i < _f.size(); i++) {
        delta_phi_x += _alpha * _f[i] * _u[i];
    }
    return delta_phi_x;
}

void FractionalPacking::update(int demand_index, const Flow &old_fxi, const Flow &flow_x_i){
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    //double low = 1/pow(_alpha, 3);
    //double low = 1/20.0/_alpha/_alpha/(_rou+_alpha);
    vector<double> ax(change_edges.size() + 1, 0);
    double ax_m = _u.back();
    vector<double> ax_star(change_edges.size() + 1, 0);
    double ax_star_m = _u.back();
    int index = 0;
    for (const auto &i: change_edges) {
        ax[index] = _u[i];
        ax_star[index++] = _u[i] + bw_change[i] * inverse_capacity[i];
        ax_star_m += bw_change[i] * beta[i];
    }
    ax[ax.size() - 1] = ax_m;
    ax_star[ax.size() - 1] = ax_star_m;
    double theta = compute_theta_newton_raphson(ax, ax_star, 0.01, 0, 1.0);
    Flow target_fxi = update_flow(old_fxi, flow_x_i, theta);
    solution.rm_flow(demand_index, bw_change, change_edges);
    solution.add_flow(demand_index, target_fxi, bw_change, change_edges);
    double old_potential = _potential;
    double new_potential = update_potential_function();
    if (new_potential < old_potential) {
        non_update=0;
        improve_flag=true;
        update_mab(demand_index, 1 - new_potential / old_potential);
        update_count++;
    } else {
        non_update++;
        if(improve_flag){
            for(int i=0; i<improve.size();i++){
                improve[i]=true;
            }
            improve_flag=false;
        }
        improve[demand_index]=false;
        update_mab(demand_index, 0);
        solution.rm_flow(demand_index, bw_change, change_edges);
        solution.add_flow(demand_index, old_fxi, bw_change, change_edges);
        update_potential_function();
    }
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        update_time += time_span.count() / 1000;
    }

}

void FractionalPacking::iteration() {
    /*
     * at each iteraion, compute the new cost, then choose an x_i to optimize, then update x..
        according to section 2.2 and section 3.1
     */

    //NOTE update_potential_function has to be called before iteration

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    //equition(6) delta phi x dot x

    /*
    double tmp1 = _alpha * _f.back();
    double tmp2 = delta_phi_x / 9.0;
    for (int i = 0; i < cost_map.size(); i++) {
        //delta_x\PHI(x); m_th row of first term; second term
        dual_cost->operator[](graph.arcFromId(i)) = int((_alpha * _f[i] * inverse_capacity[i] +
                                                         tmp1 * beta[i] + tmp2 * beta[i]));
    }
     */
    double maxcost=0;
    int SCALE=100;
    for (int i = 0; i < cost_map.size(); i++) {
        //delta_x\PHI(x); m_th row of first term; second term
        double c = _y[i]*inverse_capacity[i] + _y.back()*beta[i];
        maxcost=c>maxcost?c:maxcost;
        dual_cost->operator[](graph.arcFromId(i)) = int(_y[i]*inverse_capacity[i] + _y.back()*beta[i]);
        //according to the paper, "fast approximation algorithms for multicommodity problem. The capacity constraint should be \lambda * u, since Ax<lambda* u, constraint u cannot ensure the second inequality."
        relax_cap->operator[](graph.arcFromId(i)) = int(_rou * capacity_map[i]*SCALE);
    }
    // let x = x_1 * x_2 ... * x_k, choose x_i in round-robin order, update x_i
    
    //DEV NOTE: add a for loop here to compute min_cost_flow in a batch
    int demand_index = draw_demand_index();
    Demand demand_i = demands[demand_index];
    Flow flow_x_i = min_cost_flow(demand_i.src, demand_i.dst, demand_i.val,SCALE, dual_cost, relax_cap);
    Flow old_fxi = solution.rm_flow(demand_index, bw_change, change_edges);
    solution.add_flow(demand_index, flow_x_i, bw_change, change_edges);
    // if flow_x_i and old_fxi are the same, then no update
    bool flow_change = false;
    if (flow_x_i.size() != old_fxi.size()) {
        flow_change = true;
    } else {
        for (const auto &kv: old_fxi) {
            if (flow_x_i.count(kv.first) == 0) {
                flow_change = true;
                break;
            } else if (abs(flow_x_i[kv.first] - old_fxi[kv.first]) > 1e-5) {
                flow_change = true;
                break;
            }
        }
    }
    if (flow_change) {
        //DEV NOTE: accumulate changes in parallel but dont update the solution here
        update(demand_index, old_fxi, flow_x_i);
        if (time_debug) {
            high_resolution_clock::time_point t3 = high_resolution_clock::now();
            duration<double, std::micro> time_span = t3 - t2;
            iteration_time += time_span.count() / 1000;
        }
    } else {
        non_update++;
        update_mab(demand_index, 0);
        //NOTE lazy update on improve_flag, if the update on demand_i improves the potential, then there should be some
        if(improve_flag){
            for(int i=0; i<improve.size();i++){
                improve[i]=true;
            }
            improve_flag=false;
        }
        improve[demand_index]=false;
        if (time_debug) {
            high_resolution_clock::time_point t3 = high_resolution_clock::now();
            duration<double, std::micro> time_span = t3 - t2;
            iteration_time += time_span.count() / 1000;
        }
    }
}

int FractionalPacking::iteration_all() {
    /*
     * at each iteraion, compute the new cost, then choose an x_i to optimize, then update x..
        according to section 2.2 and section 3.1
     *
     * Return 0 if succeed
     * return 1 for early termination
     * return 2 for no improvement
     */

    //NOTE update_potential_function has to be called before iteration

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    bool res = false;
    for(int i=0; i<improve.size();i++){
        res |= improve[i];
    }
    if(!res){
        cout<<"no improvement, fail"<<endl;
        return 2;
    }
    //equition(6) delta phi x dot x
    double cost = 0;
    double solution_dual_cost = 0;
    int SCALE=100;
    for (int i = 0; i < cost_map.size(); i++) {
        //delta_x\PHI(x); m_th row of first term; second term
        double c = int(_y[i]*inverse_capacity[i] + _y.back()* beta[i]);
        dual_cost->operator[](graph.arcFromId(i)) = c;
        solution_dual_cost += solution.used_bw[i] * c;

        //+1 to make sure that relaxed capacity >=lamda* u
        relax_cap->operator[](graph.arcFromId(i)) = int(_rou * capacity_map[i]*SCALE);
    }
    FlowSolution new_sol(cost_map.size());
    for (int i = 0; i < demands.size(); i++) {
        if(i==0) {
            Flow flow = min_cost_flow(demands[i].src, demands[i].dst, demands[i].val, SCALE, dual_cost, relax_cap);
            new_sol.add_flow(i, flow);
            for (const auto &kv:flow) {
                cost += dual_cost->operator[](graph.arcFromId(kv.first)) * kv.second;
            }
        }else{
            Flow flow = min_cost_flow(demands[i].src, demands[i].dst, demands[i].val, SCALE, 0, 0);
            new_sol.add_flow(i, flow);
            for (const auto &kv:flow) {
                cost += dual_cost->operator[](graph.arcFromId(kv.first)) * kv.second;
            }
        }
    }
    //compute sum_y
    double sum_y = 0;
    for (int i = 0; i < _u.size(); i++) {
        sum_y+=_y[i];
    }

    /*
     * According to the paper lambda* sum_y >=solution cost of dual edges
     *
     */
    if (cost / sum_y > 1.00) {
        cout << "cost(f_start): " << cost << " sum_y" << sum_y * _rou
             << " solution cost: " << solution_dual_cost
             << " rou: " << _rou << endl;
        cout << cost / sum_y << endl;
        return 1;
    }
    double ratio = (sum_y * _rou / solution_dual_cost - 1) / (solution_dual_cost / cost - 1);
    if(alpha_flag) {
        //942605 before adaptive alpha 198secs
        cout << "ratio" << ratio << endl;
        if (ratio < 0.5 && _adaptive_alpha > 0) {
            _alpha = _alpha * 0.618;
            _adaptive_alpha--;
        } else{
            _alpha = theoritical_alpha;
            alpha_flag = false;
        }
        compute_potential_function(false);
    }
    assert(sum_y*_rou >= solution_dual_cost);
    assert(solution_dual_cost >=cost);

    //update
    /*
    vector<double>ax = _u;
    vector<double>ax_star(_u.size(), 0);
    double ax_star_m=0;
    for(int i=0; i<cost_map.size(); i++){
        ax_star[i] = new_sol.used_bw[i]*inverse_capacity[i];
        ax_star_m += new_sol.used_bw[i]*beta[i];
    }
    ax_star[ax_star.size()-1] = ax_star_m;

    double old_potential = _potential;
    FlowSolution old_solution = solution;
    double theta = compute_theta_newton_raphson(ax, ax_star, 0.01, 0, 1.0);
    solution.update(new_sol, theta);
    double new_potential = compute_potential_function(true);
    if(new_potential>old_potential){
        solution = old_solution;
        compute_potential_function(true);
    }
     */
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        iteration_all_time += time_span.count() / 1000;
    }
    return 0;
}

Flow FractionalPacking::update_flow(const Flow &oldx, const Flow &newx, double theta) {
    Flow res;
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    if (theta < 1) {
        for (const auto &kv: oldx) {
            res[kv.first] = (1 - theta) * kv.second;
        }
    }
    if (theta > 0) {
        for (const auto &kv: newx) {
            if (res.count(kv.first) > 0) {
                res[kv.first] += theta * kv.second;
            } else {
                res[kv.first] = theta * kv.second;
            }
        }
    }
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        update_flow_time += time_span.count() / 1000;
    }
    return res;
}

double FractionalPacking::compute_theta_newton_raphson(vector<double> &ax, vector<double> &ax_star, double theta0 =
0.5,
                                                       double mintheta = 0, double maxtheta = 1) {
    /*
     * Newton's method:
     * x1 = x0-f(x0)/f'(x0)
     */
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    double theta = theta0;
    double new_theta = update_theta(theta, ax, ax_star);
    try {
        while (new_theta < 2 && new_theta > -1 && abs(new_theta - theta) > 1e-2) {
            theta = new_theta;
            new_theta = update_theta(theta, ax, ax_star);
        }
    } catch (Exception e) {
        cout << e.what();
        new_theta = mintheta;
    }
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        new_ton_time += time_span.count() / 1000;
    }
    if (new_theta < mintheta) {
        return mintheta;
    } else if (new_theta > maxtheta) {
        return maxtheta;
    } else {
        return new_theta;
    }
}

inline double FractionalPacking::update_theta(double theta, const vector<double> &u_0, const vector<double> &u_1) {
    /*
     * delta(\phi)/delta(\theta) = sum(\alpha * f(u_i)*a_i*(x1-x0))= sum(\alpha* f*(u_i_theta0)*(u_1 -u_0),
     * where x0 is previous solution and x1 is new solution and x = (1-theta)* x0 + theta*x1
     * u_i = (1-theta)*u_0 + theta * u_1
     * we f
     */
    double f_theta0 = 0;
    double fprime_theta0 = 0;
    vector<double> t(u_0.size(), 0);
    vector<double> t1(u_1.size(), 0);
    for (int i = 0; i < u_0.size(); i++) {
        t[i] = exp(_alpha * ((1 - theta) * u_0[i] + theta * u_1[i] - 1)) * _alpha * (u_1[i] - u_0[i]);
        f_theta0 += t[i];
        t1[i] = t[i] * _alpha * (u_1[i] - u_0[i]);
        fprime_theta0 += t1[i];
    }
    return theta - f_theta0 / fprime_theta0;
}

double FractionalPacking::get_cost() {
    double c = 0;
    for (int i = 0; i < cost_map.size(); i++) {
        c += cost_map[i] * solution.used_bw[i];
    }
    return c;
}

int FractionalPacking::draw_demand_index() {
    if (mab_flag) {
        int res = _i;
        if (_i == demands.size() - 1) {
            mab_flag = false;
        }
        _i = (++_i) % demands.size();
        return res;
    } else {
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        int res = 0;
        for (_i = (++_i) % demands.size();; _i = (++_i) % demands.size()) {
            if (mab_index[_i] > mab_average) {
                if (time_debug) {
                    high_resolution_clock::time_point t3 = high_resolution_clock::now();
                    duration<double, std::micro> time_span = t3 - t2;
                    draw_index_time += time_span.count() / 1000;
                }
                return _i;
            } else if (rand() % 4 == 0) {
                if (time_debug) {
                    high_resolution_clock::time_point t3 = high_resolution_clock::now();
                    duration<double, std::micro> time_span = t3 - t2;
                    draw_index_time += time_span.count() / 1000;
                }
                return _i;
            }
        }
    }
}

void FractionalPacking::update_mab(int index, double r) {
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    //mab_reward[index] = (mab_reward[index]*mab_times[index] + r)/(++mab_times[index]);
    //mab_index[index] = mab_reward[index] + sqrt(3/mab_times[index]);
    double theta = 0;
    if(r>0){
        mab_times[index]++;
    }
    mab_average = mab_average + ((1-theta)*r - theta*mab_index[index]) / mab_index.size();
    mab_index[index] = theta*r + (1-theta)*mab_index[index];
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        draw_index_time += time_span.count() / 1000;
    }
}

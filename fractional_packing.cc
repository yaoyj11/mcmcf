//
// Created by yaoyj11 on 4/5/18.
//

#include<fstream>
#include <iostream>
#include "fractional_packing.h"
using namespace std::chrono;

FractionalPacking::FractionalPacking() :budget(0), _i(0), min_cost_count(0), capacity_map(), cost_map(), demands()
        , edges(), beta(), solution(0){
    _m = 0;
    _epsilon = 2.0;
    _alpha = 0;
    _potential = 0;
    _rou=0;
    update_count=0;
}

FractionalPacking::FractionalPacking(string filename) :budget(0), _i(0), min_cost_count(0), capacity_map(),
                                                       cost_map(), demands(), edges(), beta(), solution(0){
    update_count=0;
    _m = 0;
    _epsilon = 2.0;
    _alpha = 0;
    _potential = 0;
    _rou=0;
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
        cost_map = vector<double>(n, 0);
        for (int i = 0; i < n; i++) {
            int index, s, t, cap, cost;
            is >> index >> s >> t >> cap >> cost;
            ListDigraph::Arc arc = graph.addArc(graph.nodeFromId(s), graph.nodeFromId(t));
            capacity_map[i] = cap;
            cost_map[i] = cost;
            edges.push_back(Edge(s,t,cap,cost));
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
        cost_map = vector<double>(n, 0);
        edges = vector<Edge>();
        solution.set_arc_num(n);
        for (int i = 0; i < n; i++) {
            int index, s, t, cap, cost;
            is >> index >> s >> t >> cap >> cost;
            ListDigraph::Arc arc = graph.addArc(graph.nodeFromId(s), graph.nodeFromId(t));
            capacity_map[i] = cap;
            cost_map[i] = cost;
            edges.push_back(Edge(s,t,cap,cost));
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

Flow FractionalPacking::min_cost_flow(int src, int dst, int d, vector<double> cost = vector<double>(),
                                      vector<int> cap = vector<int>()) {
    min_cost_count++;
    //TODO: Optimize min_cost_flow
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    NetworkSimplex<ListDigraph> ns(graph);
    ListDigraph::ArcMap<int> c(graph), u(graph);
    double scale = 10.0;
    for (int i = 0; i < cost.size(); i++) {
        c[graph.arcFromId(i)] = int(cost[i] * scale);
    }
    ns.costMap(c);
    if (!cap.empty()) {
        for (int i = 0; i < cap.size(); i++) {
            u[graph.arcFromId(i)] = cap[i];
        }
        ns.upperMap(u);
    }
    ns.stSupply(graph.nodeFromId(src), graph.nodeFromId(dst), d);
    ns.run();
    Flow map;
    ListDigraph::ArcMap<int> flow_map(graph);
    ns.flowMap(flow_map);
    for (int i = 0; i < graph.maxArcId(); i++) {
        if (flow_map[graph.arcFromId(i)] > 0) {
            map[i] = flow_map[graph.arcFromId(i)];
        }
    }
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        min_cost_time+=time_span.count()/1000;
    }
    return map;
}

void FractionalPacking::set_buget(double b) {
    budget = b;
    beta = vector<double>(edges.size(), 0);
    for(int i = 0; i< edges.size(); i++){
        beta[i] = cost_map[i]/budget;
    }
}

FlowSolution FractionalPacking::fractional_packing(double epsilon=0.05, bool restart) {
    _m = edges.size() + 1;
    _u = vector<double>(_m,0);
    _f = vector<double>(_m,0);
    change = vector<bool>(edges.size(), false);
    bw_change = vector<double>(edges.size(),0);
    _epsilon = edges.size() -1;
    if(time_debug){

        init_flow_time = 0;

        min_cost_time=0;

        new_ton_time=0;

        potential_time=0;

        update_flow_time=0;

        iteration_time=0;
    }
    if(restart){
        compute_init_flow();
    }
    while(_epsilon>=epsilon){
        cout<<current_date_time()<< " epsilon: "<<_epsilon<<endl;
        compute_potential_function();
        while(_potential > 3 * _m) {
            //while(_rou>1+_epsilon){
            iteration();
           _i = (_i+1)%demands.size();
        //cout<<potential<<endl;
        }
        if(_epsilon == epsilon){
            break;
        }
        _epsilon /= 2.0;
        if(_epsilon<epsilon){
            _epsilon = epsilon;
        }
    }
    return solution;
}

void FractionalPacking::compute_init_flow() {
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    solution.set_arc_num(edges.size());
    for(int i=0;i<demands.size();i++){
        Demand d = demands[i];
        Flow f = min_cost_flow(d.src, d.dst, d.val, cost_map, capacity_map);
        if(f.empty()){
            throw std::runtime_error("No feasible solution for single commodity");
        }
        solution.add_flow(i, f, change, bw_change);
    }
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        init_flow_time+=time_span.count()/1000;
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

double FractionalPacking::compute_potential_function() {
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    _alpha = 1/_epsilon*log(3*_m);
    for(int i = 0; i< edges.size(); i++){
        if(change[i]) {
            _u[i] = solution.used_bw[i] / capacity_map[i];
            _u[_u.size() - 1] += bw_change[i] * beta[i];
            change[i] = false;
            bw_change[i] = 0;
        }
    }
    _rou=0;
    for(int i=0;i<_u.size(); i++){
        _rou = _rou>_u[i]?_rou:_u[i];
    }
    _potential = 0;
    for(int i = 0; i<_u.size(); i++){
        _f[i] = exp(_alpha * (_u[i] - 1.0));
        _potential += _f[i];
    }
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        potential_time+=time_span.count()/1000;
    }
    return _potential;
}

double FractionalPacking::update_potential_function() {
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    _alpha = 1/_epsilon*log(3*_m);
    for(int i = 0; i< edges.size(); i++){
        if(change[i]) {
            _u[i] = solution.used_bw[i] / capacity_map[i];
            _u[_u.size() - 1] += bw_change[i] * beta[i];
        }
    }
    _rou=0;
    for(int i=0;i<_u.size(); i++){
        _rou = _rou>_u[i]?_rou:_u[i];
    }
    _potential = 0;
    for(int i = 0; i<change.size(); i++){
        if(change[i]) {
            _f[i] = exp(_alpha * (_u[i] - 1.0));
            change[i]=false;
            bw_change[i]=0;
        }
        _potential += _f[i];
    }
    _f[_f.size()-1]=exp(_alpha*(_u.back() - 1.0));
    _potential += _f.back();
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        potential_time+=time_span.count()/1000;
    }
    return _potential;
}

void FractionalPacking::iteration() {
    /*
     * at each iteraion, compute the new cost, then choose an x_i to optimize, then update x..
        according to section 2.2 and section 3.1
     */

    //NOTE update_potential_function has to be called before iteration

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    //equition(6) delta phi x dot x
    double delta_phi_x = 0;

    //compute new cost in (4)

    for (int i = 0; i< _f.size(); i++){
        delta_phi_x += _alpha* _f[i]* _u[i];
    }

    vector<double> dual_cost(edges.size(), 0);
    for(int i = 0; i<edges.size(); i++){
        //delta_x\PHI(x); m_th row of first term; second term
        dual_cost[i] = _alpha* _f[i]/capacity_map[i] + _alpha*_f[_f.size()-1] * beta[i] + delta_phi_x * beta[i]
                                                                                        /9.0/_alpha;
    }
    // let x = x_1 * x_2 ... * x_k, choose x_i in round-robin order, update x_i
    Demand demand_i = demands[_i];
    Flow flow_x_i = min_cost_flow(demand_i.src, demand_i.dst, demand_i.val,dual_cost, capacity_map);
    int m=solution.flows.size();
    Flow old_fxi = solution.rm_flow(_i, change, bw_change);
    solution.add_flow(_i, flow_x_i, change, bw_change);
    //double low = 1/pow(_alpha, 3);
    double low = 1/20.0/_alpha/_alpha/(_rou+_alpha);
    // if flow_x_i and old_fxi are the same, then no update
    //TODO: if no update, there is no need to compute potential_function again
    bool flow_change = false;
    if(flow_x_i.size()!= old_fxi.size()){
        flow_change = true;
    }else{
        for(const auto&kv: old_fxi){
            if(flow_x_i.count(kv.first)==0){
                flow_change = true;
                break;
            }else if(abs(flow_x_i[kv.first] - old_fxi[kv.first])>1e-5){
                flow_change = true;
                break;
            }
        }
    }
    if(flow_change){
        double theta = compute_theta_newton_raphson(low, 0, 1.0);
        Flow target_fxi = update_flow(old_fxi, flow_x_i, theta);
        solution.rm_flow(_i, change, bw_change);
        solution.add_flow(_i, target_fxi, change, bw_change);
        double old_potential = _potential;
        double new_potential = update_potential_function();
        if(new_potential<old_potential){
            update_count++;
            if (time_debug) {
                high_resolution_clock::time_point t3 = high_resolution_clock::now();
                duration<double, std::micro> time_span = t3 - t2;
                iteration_time+=time_span.count()/1000;
            }
            return;
        }else{
            solution.rm_flow(_i, change, bw_change);
            solution.add_flow(_i, old_fxi, change, bw_change);
            update_potential_function();
            if (time_debug) {
                high_resolution_clock::time_point t3 = high_resolution_clock::now();
                duration<double, std::micro> time_span = t3 - t2;
                iteration_time+=time_span.count()/1000;
            }
        }
    }
}

Flow FractionalPacking::update_flow(Flow oldx, Flow newx, double theta) {
    Flow res;
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    if(theta<1) {
        for (const auto &kv: oldx) {
            res[kv.first] = (1 - theta) * kv.second;
        }
    }
    if(theta>0) {
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
        update_flow_time+=time_span.count()/1000;
    }
    return res;
}

double FractionalPacking::compute_theta_newton_raphson(double theta0 = 0.5,
                                                       double mintheta=0, double maxtheta=1) {
    /*
     * Newton's method:
     * x1 = x0-f(x0)/f'(x0)
     */
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    vector<double> ax;
    double ax_m = _u.back();
    vector<double>ax_star;
    double ax_star_m = _u.back();
    for (int i = 0; i< edges.size(); i++){
        if(change[i]) {
            ax.push_back( _u[i]);
            ax_star.push_back(_u[i] + bw_change[i] / capacity_map[i]);
            ax_star_m += bw_change[i] * beta[i];
        }
    }
    ax.push_back(ax_m);
    ax_star.push_back(ax_star_m);
    double theta = theta0;
    double new_theta = update_theta(theta, ax, ax_star);
    try{
        while(new_theta<2&&new_theta>-1 &&abs(new_theta - theta)>1e-3){
            theta = new_theta;
            new_theta = update_theta(theta, ax, ax_star);
        }
    }catch (Exception e){
        cout<<e.what();
        new_theta = mintheta;
    }
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        new_ton_time+=time_span.count()/1000;
    }
    if(new_theta < mintheta){
        return mintheta;
    }else if(new_theta>maxtheta){
        return maxtheta;
    }else{
        return new_theta;
    }
}

double FractionalPacking::update_theta(double theta, vector<double> u_0, vector<double> u_1) {
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
    double p=0;
    for(int i = 0; i< u_0.size(); i++){
        t[i] = exp(_alpha* ((1-theta)*u_0[i] + theta* u_1[i] - 1)) * _alpha * (u_1[i]-u_0[i]);
        p+=exp(_alpha*((1-theta)*u_0[i]+theta*u_1[i]-1));
        f_theta0 += t[i];
        t1[i] = t[i]*_alpha*(u_1[i]-u_0[i]);
        fprime_theta0 += t1[i];
    }
    return theta - f_theta0/fprime_theta0;
}

double FractionalPacking::get_cost() {
    double c =0;
    for(int i = 0; i<edges.size(); i++){
        c += cost_map[i] * solution.used_bw[i];
    }
    return c;
}


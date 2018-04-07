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
}

FractionalPacking::FractionalPacking(string filename) :budget(0), _i(0), min_cost_count(0), capacity_map(),
                                                       cost_map(), demands(), edges(), beta(), solution(0){
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
        std::cout << "min cost flow " << time_span.count() / 1000 << " milli seconds" << endl;
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
    _delta_phi = vector<double>(edges.size(),0);
    _epsilon = edges.size() -1;
    time_debug=false;
    if(restart){
        compute_init_flow();
    }
    while(_epsilon>=epsilon){
        time_debug=false;
        cout<<current_date_time()<< " epsilon: "<<_epsilon<<endl;
        double potential = update_potential_function();
        while(potential > 3 * _m) {
            //while(_rou>1+_epsilon){
            iteration();
            time_debug=false;
           _i = (_i+1)%demands.size();
           potential = update_potential_function();
           cout<<potential<<endl;
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
    solution.set_arc_num(edges.size());
    for(int i=0;i<demands.size();i++){
        Demand d = demands[i];
        Flow f = min_cost_flow(d.src, d.dst, d.val, cost_map, capacity_map);
        if(f.empty()){
            throw std::runtime_error("No feasible solution for single commodity");
        }
        solution.add_flow(i, f);
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

double FractionalPacking::compute_potential_function(FlowSolution sol) {
    vector<double> u(_m,0);
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    double cost = 0;
    for(int i = 0; i< edges.size(); i++){
        double bw = sol.used_bw[i];
        u[i] = bw/capacity_map[i];
        cost+=bw*cost_map[i];
    }
    u[u.size()-1]=cost/budget;
    double potential = 0;
    for(int i = 0; i<u.size(); i++){
        if(_alpha*(u[i] - 1) < 20) {
            potential += exp(_alpha * (u[i] - 1.0));
        }else{
            potential+=10000;
        }
    }
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        std::cout << "compute_potential_function " << time_span.count() / 1000 << " milli seconds" << endl;
    }
    return potential;
}

double FractionalPacking::update_potential_function() {
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    _alpha = 1/_epsilon*log(3*_m);
    double cost = 0;
    for(int i = 0; i< edges.size(); i++){
        double bw = solution.used_bw[i];
        _u[i] = bw/capacity_map[i];
        cost+=bw*cost_map[i];
    }
    _u[_u.size()-1]=cost/budget;
    _rou=0;
    for(int i=0;i<_u.size(); i++){
        _rou = _rou>_u[i]?_rou:_u[i];
    }
    _potential = 0;
    for(int i = 0; i<_u.size(); i++){
        _f[i] = exp(_alpha * (_u[i] - 1.0));
        _potential += _f[i];
    }
    //compute delta_phi
    for(int i = 0; i< edges.size(); i++){
        _delta_phi[i] = _alpha * _f[i] / capacity_map[i] + _alpha * _f[_f.size()-1]*beta[i];
    }
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        std::cout << "update_potential_function " << time_span.count() / 1000 << " milli seconds" << endl;
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
    vector<double> x0 = solution.used_bw;
    Flow old_fxi = solution.rm_flow(_i);
    //TODO: problem using demand_i as key
    solution.add_flow(_i, flow_x_i);
    vector<double> x1= solution.used_bw;
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
            }else if(abs(flow_x_i[kv.first] - old_fxi[kv.first])>1e-6){
                //cout<<"flow diff "<<kv.first<<" "<<flow_x_i[kv.first]<<" "<<old_fxi[kv.first]<<endl;
                flow_change = true;
                break;
            }
        }
    }
    if(flow_change){
        double theta = compute_theta_newton_raphson(x0, x1, low, low, 1.0);

        Flow target_fxi = update_flow(old_fxi, flow_x_i, theta);
        solution.rm_flow(_i);
        solution.add_flow(_i, target_fxi);
        double new_potential = compute_potential_function(solution);
        if(new_potential<_potential){
            if (time_debug) {
                high_resolution_clock::time_point t3 = high_resolution_clock::now();
                duration<double, std::micro> time_span = t3 - t2;
                std::cout << "iteration " << time_span.count() / 1000 << " milli seconds" << endl;
            }
            return;
        }else{
            solution.rm_flow(_i);
            solution.add_flow(_i, old_fxi);
            if (time_debug) {
                high_resolution_clock::time_point t3 = high_resolution_clock::now();
                duration<double, std::micro> time_span = t3 - t2;
                std::cout << "iteration " << time_span.count() / 1000 << " milli seconds" << endl;
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
        std::cout << "update_flow " << time_span.count() / 1000 << " milli seconds" << endl;
    }
    return res;
}

double FractionalPacking::compute_theta_newton_raphson(vector<double> x0, vector<double> x1, double theta0 = 0.5,
                                                       double mintheta=0, double maxtheta=1) {
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    vector<double> ax(_m, 0);
    double ax_m = 0;
    vector<double>ax_star(_m, 0);
    double ax_star_m = 0;
    for (int i = 0; i< edges.size(); i++){
        ax[i] = x0[i]/capacity_map[i];
        ax_m += x0[i] * cost_map[i]/budget;
        ax_star[i] = x1[i]/capacity_map[i];
        ax_star_m += x1[i]*cost_map[i]/budget;
    }
    ax[_m - 1] = ax_m;
    ax_star[_m-1] = ax_star_m;
    double theta = theta0;
    double new_theta = update_theta(theta, ax, ax_star);
    try{
        while(abs(new_theta - theta)>1e-2){
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
        std::cout << "compute_theta_newton_raphson " << time_span.count() / 1000 << " milli seconds" << endl;
    }
    if(new_theta < mintheta){
        return mintheta;
    }else if(new_theta>maxtheta){
        return maxtheta;
    }else{
        return new_theta;
    }
}

double FractionalPacking::update_theta(double theta, vector<double> ax, vector<double> ax_star) {
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    double f_theta0 = 0;
    double fprime_theta0 = 0;
    vector<double> t(ax.size(), 0);
    vector<double> t1(ax_star.size(), 0);
    for(int i = 0; i< ax.size(); i++){
        t[i] = exp(_alpha* (ax[i] + theta* (ax_star[i] - ax[i]) - 1))*_alpha*(ax_star[i]-ax[i]);
        f_theta0 += t[i];
        t1[i] = t[i]*_alpha*(ax_star[i]-ax[i]);
        fprime_theta0 += t1[i];
    }
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        std::cout << "update_theta " << time_span.count() / 1000 << " milli seconds" << endl;
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


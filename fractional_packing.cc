//
// Created by yaoyj11 on 4/5/18.
//

#include<fstream>
#include <iostream>
#include "fractional_packing.h"
using namespace std::chrono;

FractionalPacking::FractionalPacking() :budget(0), _i(0), min_cost_count(0), capacity_map(), cost_map(), demands()
        , beta(), solution(0){
    _m = 0;
    _epsilon = 2.0;
    _alpha = 0;
    _potential = 0;
    _rou=0;
    update_count=0;
    network_simplex=0;
    dual_cost=0;
    delta_phi_x=0;
}

FractionalPacking::FractionalPacking(string filename) :budget(0), _i(0), min_cost_count(0), capacity_map(),
                                                       cost_map(), demands(), beta(), solution(0){
    network_simplex=0;
    dual_cost=0;
    update_count=0;
    _m = 0;
    _epsilon = 2.0;
    _alpha = 0;
    _potential = 0;
    _rou=0;
    delta_phi_x=0;
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
            inverse_capacity[i] = 1.0/cap;
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
    if(network_simplex!=0){
        delete network_simplex;
        network_simplex=0;
    }
    if(dual_cost!=0){
        delete dual_cost;
        dual_cost=0;
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
        inverse_capacity = vector<double>(n,0);
        cost_map = vector<double>(n, 0);
        solution.set_arc_num(n);
        for (int i = 0; i < n; i++) {
            int index, s, t, cap, cost;
            is >> index >> s >> t >> cap >> cost;
            ListDigraph::Arc arc = graph.addArc(graph.nodeFromId(s), graph.nodeFromId(t));
            capacity_map[i] = cap;
            inverse_capacity[i] = 1.0/cap;
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

FlowSolution FractionalPacking::fractional_packing(double epsilon=0.05, bool restart) {
    _m = cost_map.size() + 1;
    _u = vector<double>(_m,0);
    _f = vector<double>(_m,0);
    mab_reward = vector<double>(demands.size(),0);
    mab_index = vector<double>(demands.size(), 0);
    mab_times = vector<int>(demands.size(), 0);
    mab_average=0;
    mab_flag=true;
    bw_change = vector<double>(cost_map.size(),0);
    _epsilon = cost_map.size() -1;
    if(time_debug){

        init_flow_time = 0;

        min_cost_time=0;

        new_ton_time=0;

        potential_time=0;

        update_flow_time=0;

        iteration_time=0;

        draw_index_time = 0;
    }
    if(restart){
        network_simplex = new NetworkSimplex<ListDigraph>(graph);
        dual_cost = new ListDigraph::ArcMap<int>(graph);
        compute_init_flow();
        compute_potential_function(true);
    }
    while(_epsilon>=epsilon){
        cout<<current_date_time()<< " epsilon: "<<_epsilon<<endl;
        compute_potential_function();
        while(_potential > 3 * _m) {
            //while(_rou>1+_epsilon){
            iteration();
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

Flow FractionalPacking::min_cost_flow(int src, int dst, int d, const vector<double> &cost,
                                      vector<int> cap) {
    min_cost_count++;
    //TODO: Optimize min_cost_flow
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    ListDigraph::ArcMap<int> c(graph);
    double scale = 10.0;
    for (int i = 0; i < cost.size(); i++) {
        c[graph.arcFromId(i)] = int(cost[i] * scale);
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
        min_cost_time+=time_span.count()/1000;
    }
    return map;
}

Flow FractionalPacking::min_cost_flow(int src, int dst, int d, ListDigraph::ArcMap<int> *c) {
    min_cost_count++;
    //TODO: Optimize min_cost_flow
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    if(c!=0) {
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
        min_cost_time+=time_span.count()/1000;
    }
    return map;
}

void FractionalPacking::set_buget(double b) {
    budget = b;
    beta = vector<double>(cost_map.size(), 0);
    for(int i = 0; i< cost_map.size(); i++){
        beta[i] = cost_map[i]/budget;
    }
}


void FractionalPacking::compute_init_flow() {
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    solution.set_arc_num(cost_map.size());
    for(int i=0;i<demands.size();i++){
        Demand d = demands[i];
        Flow f = min_cost_flow(d.src, d.dst, d.val, cost_map, capacity_map);
        if(f.empty()){
            throw std::runtime_error("No feasible solution for single commodity");
        }
        solution.add_flow(i, f, bw_change, change_edges);
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

double FractionalPacking::compute_potential_function(bool recompute_u) {
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    _alpha = 1/_epsilon*log(3*_m);
    if(recompute_u){
        _u[_u.size()-1]=0;
        for(int i=0; i<cost_map.size(); i++) {
            _u[i] = solution.used_bw[i] * inverse_capacity[i];
            _u[_u.size() - 1] += solution.used_bw[i] * beta[i];
            bw_change[i] = 0;
        }
    }
    change_edges.clear();
    /*
    _rou=0;
    for(int i=0;i<_u.size(); i++){
        _rou = _rou>_u[i]?_rou:_u[i];
    }
     */
    _potential = 0;
    delta_phi_x=0;
    for(int i = 0; i<_u.size(); i++){
        _f[i] = exp(_alpha * (_u[i] - 1.0));
        _potential += _f[i];
        delta_phi_x += _alpha* _f[i]* _u[i];
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
    for(const auto&i:change_edges){
        _u[i] = solution.used_bw[i] *inverse_capacity[i];
        _u[_u.size() - 1] += bw_change[i] * beta[i];
        bw_change[i]=0;
    }
    /*
    _rou=0;
    for(int i=0;i<_u.size(); i++){
        _rou = _rou>_u[i]?_rou:_u[i];
    }*/
    _potential -= _f.back();
    delta_phi_x -= _alpha*_f.back()*_u.back();
    for(const auto&i:change_edges){
        _potential -=_f[i];
        delta_phi_x-=_alpha*_f[i]*_u[i];
        _f[i] = exp(_alpha * (_u[i] - 1.0));
        _potential += _f[i];
        delta_phi_x +=_alpha*_f[i]*_u[i];
    }
    _f[_f.size()-1]=exp(_alpha*(_u.back() - 1.0));
    _potential += _f.back();
    delta_phi_x += _alpha*_f.back()*_u.back();
    change_edges.clear();
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        potential_time+=time_span.count()/1000;
    }
    return _potential;
}

double FractionalPacking::compute_delta_phi_x(){
    double delta_phi_x = 0;
    //compute new cost in (4)
    for (int i = 0; i< _f.size(); i++){
        delta_phi_x += _alpha* _f[i]* _u[i];
    }
    return delta_phi_x;
}

void FractionalPacking::iteration() {
    /*
     * at each iteraion, compute the new cost, then choose an x_i to optimize, then update x..
        according to section 2.2 and section 3.1
     */

    //NOTE update_potential_function has to be called before iteration

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    //equition(6) delta phi x dot x

    double tmp1 = _alpha * _f.back();
    double tmp2 = delta_phi_x/9.0/_alpha;
    for(int i = 0; i<cost_map.size(); i++){
        //delta_x\PHI(x); m_th row of first term; second term
        dual_cost->operator[](graph.arcFromId(i)) = int((_alpha* _f[i]*inverse_capacity[i] +
                                                   tmp1 * beta[i] + tmp2 * beta[i])*10.0);
    }
    // let x = x_1 * x_2 ... * x_k, choose x_i in round-robin order, update x_i
    int demand_index = draw_demand_index();
    Demand demand_i = demands[demand_index];
    Flow flow_x_i = min_cost_flow(demand_i.src, demand_i.dst, demand_i.val,dual_cost);
    Flow old_fxi = solution.rm_flow(demand_index, bw_change, change_edges);
    solution.add_flow(demand_index, flow_x_i, bw_change, change_edges);
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
        //double low = 1/pow(_alpha, 3);
        //double low = 1/20.0/_alpha/_alpha/(_rou+_alpha);
        vector<double> ax(change_edges.size()+1, 0);
        double ax_m = _u.back();
        vector<double>ax_star(change_edges.size() +1, 0);
        double ax_star_m = _u.back();
        int index =0;
        for (const auto&i: change_edges){
            ax[index] = _u[i];
            ax_star[index++] = _u[i] + bw_change[i] *inverse_capacity[i];
            ax_star_m += bw_change[i] * beta[i];
        }
        ax[ax.size()-1]=ax_m;
        ax_star[ax.size() -1]=ax_star_m;
        double theta = compute_theta_newton_raphson(ax, ax_star,0.01, 0, 1.0);
        Flow target_fxi = update_flow(old_fxi, flow_x_i, theta);
        solution.rm_flow(demand_index, bw_change, change_edges);
        solution.add_flow(demand_index, target_fxi, bw_change, change_edges);
        double old_potential = _potential;
        double new_potential = update_potential_function();
        if(new_potential<old_potential){
            update_mab(demand_index, 1-new_potential/old_potential);
            update_count++;
            if (time_debug) {
                high_resolution_clock::time_point t3 = high_resolution_clock::now();
                duration<double, std::micro> time_span = t3 - t2;
                iteration_time+=time_span.count()/1000;
            }
            return;
        }else{
            update_mab(demand_index, 0);
            solution.rm_flow(demand_index, bw_change, change_edges);
            solution.add_flow(demand_index, old_fxi, bw_change, change_edges);
            compute_potential_function(true);
            if (time_debug) {
                high_resolution_clock::time_point t3 = high_resolution_clock::now();
                duration<double, std::micro> time_span = t3 - t2;
                iteration_time+=time_span.count()/1000;
            }
        }
    }else{
        update_mab(demand_index, 0);
    }
}

void FractionalPacking::iteration_all() {
    /*
     * at each iteraion, compute the new cost, then choose an x_i to optimize, then update x..
        according to section 2.2 and section 3.1
     */

    //NOTE update_potential_function has to be called before iteration

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    //equition(6) delta phi x dot x

    double tmp1 = _alpha * _f.back();
    double tmp2 = delta_phi_x/9.0/_alpha;
    for(int i = 0; i<cost_map.size(); i++){
        //delta_x\PHI(x); m_th row of first term; second term
        dual_cost->operator[](graph.arcFromId(i)) = int((_alpha* _f[i]*inverse_capacity[i] +
                                                         tmp1 * beta[i] + tmp2 * beta[i])*10.0);
    }
    // let x = x_1 * x_2 ... * x_k, choose x_i in round-robin order, update x_i
    FlowSolution new_solution(cost_map.size());
    int demand_index = draw_demand_index();

    for(int i=0;i<demands.size();i++){
        Demand demand_i = demands[i];
        if(i!=0) {
            new_solution.add_flow(demand_index, min_cost_flow(demand_i.src, demand_i.dst, demand_i.val),
                                  bw_change, change_edges);
        }else{
            new_solution.add_flow(demand_index, min_cost_flow(demand_i.src, demand_i.dst, demand_i.val, dual_cost),
                                  bw_change, change_edges);
        }
    }
    vector<double> ax=_u;
    vector<double> ax_star(_u.size(),0);
    double  ax_start_m = 0;
    for(int i=0; i<cost_map.size(); i++){
        ax_star[i] = new_solution.used_bw[i]*inverse_capacity[i];
        ax_start_m += new_solution.used_bw[i]*beta[i];
    }
    ax_star[ax_star.size()-1] = ax_start_m;

    double theta = compute_theta_newton_raphson(ax, ax_star,0.01, 0, 1.0);
    double old_potential = _potential;
    FlowSolution old_solution = solution;
    solution.update(new_solution, theta);
    double new_potential = compute_potential_function(true);
    if(new_potential<old_potential){
        if (time_debug) {
            high_resolution_clock::time_point t3 = high_resolution_clock::now();
            duration<double, std::micro> time_span = t3 - t2;
            iteration_time+=time_span.count()/1000;
        }
        return;
    }else{
        solution = old_solution;
        compute_potential_function(true);
        if (time_debug) {
            high_resolution_clock::time_point t3 = high_resolution_clock::now();
            duration<double, std::micro> time_span = t3 - t2;
            iteration_time+=time_span.count()/1000;
        }
    }
}

Flow FractionalPacking::update_flow(const Flow &oldx, const Flow &newx, double theta) {
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

double FractionalPacking::compute_theta_newton_raphson(vector<double> & ax, vector<double> &ax_star, double theta0 =
0.5,
                                                       double mintheta=0, double maxtheta=1) {
    /*
     * Newton's method:
     * x1 = x0-f(x0)/f'(x0)
     */
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    double theta = theta0;
    double new_theta = update_theta(theta, ax, ax_star);
    try{
        while(new_theta<2&&new_theta>-1 &&abs(new_theta - theta)>1e-2){
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
    for(int i = 0; i< u_0.size(); i++){
        t[i] = exp(_alpha* ((1-theta)*u_0[i] + theta* u_1[i] - 1)) * _alpha * (u_1[i]-u_0[i]);
        f_theta0 += t[i];
        t1[i] = t[i]*_alpha*(u_1[i]-u_0[i]);
        fprime_theta0 += t1[i];
    }
    return theta - f_theta0/fprime_theta0;
}

double FractionalPacking::get_cost() {
    double c =0;
    for(int i = 0; i<cost_map.size(); i++){
        c += cost_map[i] * solution.used_bw[i];
    }
    return c;
}

int FractionalPacking::draw_demand_index(){
    if(mab_flag) {
        int res = _i;
        if(_i==demands.size()-1){
            mab_flag= false;
        }
        _i = (++_i)%demands.size();
        return res;
    } else{
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        int res = 0;
        for(_i=(++_i)%demands.size();;_i=(++_i)%demands.size()){
            if(mab_index[_i]>mab_average){
                if (time_debug) {
                    high_resolution_clock::time_point t3 = high_resolution_clock::now();
                    duration<double, std::micro> time_span = t3 - t2;
                    draw_index_time+=time_span.count()/1000;
                }
                return _i;
            }else if(rand()%3==0){
                if (time_debug) {
                    high_resolution_clock::time_point t3 = high_resolution_clock::now();
                    duration<double, std::micro> time_span = t3 - t2;
                    draw_index_time+=time_span.count()/1000;
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
    double theta = 0.8;
    mab_average = mab_average+(r-mab_index[index])/mab_index.size();
    mab_index[index] = r;
    if (time_debug) {
        high_resolution_clock::time_point t3 = high_resolution_clock::now();
        duration<double, std::micro> time_span = t3 - t2;
        draw_index_time+=time_span.count()/1000;
    }}

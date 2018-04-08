//
// Created by yaoyj11 on 4/6/18.
//

#ifndef MCMCF_EDGE_H
#define MCMCF_EDGE_H

class Edge {
public:
    int src;
    int dst;
    int capacity;
    int cost;

    Edge(int s, int e, int cap, int c) : src(s), dst(e), capacity(cap), cost(c) {}

    Edge() {}

    Edge(const Edge &e) {
        src = e.src;
        dst = e.dst;
        capacity = e.capacity;
        cost = e.cost;
    }
};

#endif //MCMCF_EDGE_H

//
// Created by yaoyj11 on 4/5/18.
//

#ifndef MCMCF_DEMAND_H
#define MCMCF_DEMAND_H

class Demand {
public:
    int src, dst;
    int val;
    Demand(){}
    Demand(int psrc, int pdst, int pval){
        src = psrc;
        dst = pdst;
        val = pval;
    }

    bool operator ==(const Demand& d1) const{
        return src==d1.src && dst==d1.dst&&val==d1.val;
    }

    bool operator <(const Demand& d1) const{
        //This method is meaningless as defined to use Demand object as key
        return src+dst+val-d1.src-d1.dst-d1.val > 0;
    }

};

#endif //MCMCF_DEMAND_H

// ==============================================================================
//
//  RandNum.h
//  QTR
//
//  Created by Albert Lu on 8/4/18.
//  alu@tacc.utexas.edu
//
//  Last modified on 8/4/18
//
//  Note:
//
// ==============================================================================

#ifndef QTR_RANDNUM_H
#define QTR_RANDNUM_H

#include "Eigen.h"
#include "Pointers.h"

#include "./libraries/PCG/pcg-cpp-0.98/pcg_random.hpp"

namespace QTR_NS  {

    class RandNum : protected Pointers  {
        
    public:
        RandNum(class QTR *);
        ~RandNum();
        
        void     ran_mars(int seed);
        double   ran2(long &seed);
        double   mars_uniform();
        double   pcg_uniform();
        double   uniform();
        double   gaussian();
        long     randInteger(long number);
        long     getRandState();
        VectorXd   sphere();
        
    private:
        
        Error    *err;
        Log    *log;
        Parameters *parameters;
        
        long     rSeed;
        int    rank;
        
        int    save;
        int    i97;
        int    j97;
        double   *u;
        double   second;
        double   c;
        double   cd;
        double   cm;
        
        // PCG
        pcg32    rng;
    };
    
}

#endif /* QTR_RANDNUM_H */

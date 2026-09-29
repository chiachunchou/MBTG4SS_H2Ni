// ==============================================================================
//
//  RandNum.cpp
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

#include <time.h>
#include <math.h>
#include <random>

#include "Qtr.h"
#include "Error.h"
#include "Log.h"
#include "Parameters.h"
#include "RandNum.h"

using namespace QTR_NS;

const char MARS[]  = "mars";
const char PCG[]   = "pcg";
const char RAN2[]  = "ran2";

std::uniform_real_distribution<double> dist(0.0, 1.0);

/* ------------------------------------------------------------------------------- */

RandNum::RandNum(QTR *qtr) : Pointers(qtr)
{
    err = qtr->error;
    log = qtr->log;
    parameters = qtr->parameters;
    rSeed = parameters->rngSeed;
    
    if (parameters->rngType == MARS)  {
        
        ran_mars(rSeed);
    }
    else if (parameters->rngType == PCG)  {
        
        rng.seed((uint64_t)rSeed);
    }
    
    save = 0;
}
/* ------------------------------------------------------------------------------- */

RandNum::~RandNum()
{
    if (parameters->rngType == MARS)  {
        
        delete [] u;
    }
    return;
}
/* ------------------------------------------------------------------------------- */

double RandNum::uniform()
{
    if (parameters->rngType == MARS)  {
    
        return mars_uniform();
    }
    else if (parameters->rngType == PCG)  {
    
        return pcg_uniform();
    }
    else if (parameters->rngType == RAN2)  {
    
        return ran2(rSeed);
    }
    else  {
        log->log("[RandNum] Error: Invalid RNG type.\n");
        log->log("[RandNum] Please use mars, pcg, or ran2.\n");
        
        err->abort_all();
        
        return -1;
    }
}
/* ------------------------------------------------------------------------------- */

void RandNum::ran_mars(int seed)
{
    int ij, kl, i, j, k, l, ii, jj, m;
    double s, t;
    
    if (seed <= 0 || seed > 900000000)
    {
        log->log("[RandNum] Invalid seed for Marsaglia RNG (>0-900000000)\n");
        log->log("[RandNum] Generating RN seed using current time.\n");
        
        srand ((unsigned)time(NULL));
        
# ifdef QTRMPI
        MPI_Comm_rank(parameters->universe, &rank);
# else
        rank = 0;
# endif
        for (int iter = 0; iter < 10 * rank + 1; iter++)
            seed = rand() + rank * (rank + 23);
        
        seed = seed % 900000000;
        log->log("[RandNum] CPU_%d Use seed = %ld\n", rank, seed);
    }
    
    u = new double[97+1];
    
    ij = (seed - 1) / 30082;
    kl = (seed - 1) - 30082 * ij;
    i = (ij / 177) % 177 + 2;
    j = ij % 177 + 2;
    k = (kl / 169) % 178 + 1;
    l = kl % 169;
    
    for (ii = 1; ii <= 97; ii++)
    {
        s = 0.0;
        t = 0.5;
        
        for (jj = 1; jj <= 24; jj++)
        {
            m = ((i * j) % 179) * k % 179;
            i = j;
            j = k;
            k = m;
            l = (53 * l + 1) % 169;
            
            if ((l * m) % 64 >= 32) s = s + t;
            
            t = 0.5 * t;
        }
        u[ii] = s;
    }
    c = 362436.0 / 16777216.0;
    cd = 7654321.0 / 16777216.0;
    cm = 16777213.0 / 16777216.0;
    i97 = 97;
    j97 = 33;
    mars_uniform();
}
/* ------------------------------------------------------------------------------- */

double RandNum::mars_uniform()
{
    double uni = u[i97] - u[j97];
    
    if (uni < 0.0) uni += 1.0;
    
    u[i97] = uni;
    i97--;
    
    if (i97 == 0) i97 = 97;
    
    j97--;
    
    if (j97 == 0) j97 = 97;
    
    c -= cd;
    
    if (c < 0.0) c += cm;
    
    uni -= c;
    
    if (uni < 0.0) uni += 1.0;
    
    return uni;
}
/* ------------------------------------------------------------------------------- */

double RandNum::ran2(long &seed)
{
    // Long period (> 2E+18) RNG of L’Ecuyer with Bays-Durham shuffle and added safeguards.
    // Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values)
    // Call with idum a negative integer to initialize
    
    const int IM1 = 2147483563;
    const int IM2 = 2147483399;
    const int IA1 = 40014;
    const int IA2 = 40692;
    const int IQ1 = 53668;
    const int IQ2 = 52774;
    const int IR1 = 12211;
    const int IR2 = 3791;
    const int NTAB = 32;
    const int IMM1 = IM1 - 1;
    const int NDIV = 1 + IMM1 / NTAB;
    
    // RNMX should approximate the largest floating value that is less than 1
    // EPS may be calculated using double before_1 = 1 - nextafter(1.0, 0.0);
    const double EPS = 1.11E-16;
    const double RNMX = 1.0 - EPS;
    const double AM = 1.0 / IM1;
    
    static long seed2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    
    long j, k;
    double temp;
    
    if (seed <= 0)  {
        
        log->log("[RandNum] Invalid seed for RAN2 RNG (>=1)\n");
        
        log->log("[RandNum] Generating RN seed using current time.\n");
        
        srand ((unsigned)time(NULL));
        
# ifdef BSMPI
        MPI_Comm_rank(parameters->universe, &rank);
# else
        rank = 0;
# endif
        for (int iter = 0; iter < 10 * rank + 1; iter++)
            seed = rand() + rank * (rank + 23);
        
        log->log("[RandNum] CPU_%d Use seed = %ld\n", rank, seed);
    
        seed2 = seed;
        
        
        for (j = NTAB + 7; j >= 0; j--) {
            k = (seed) / IQ1;
            seed = IA1 * (seed - k * IQ1) - k * IR1;
            
            if (seed < 0) seed += IM1;
            
            if (j < NTAB) iv[j] = seed;
        }
        iy = iv[0];
    }
    k = seed / IQ1;
    seed = IA1 * (seed - k * IQ1) - k * IR1;
    
    if (seed < 0) seed += IM1;
    
    k = seed2 / IQ2;
    seed2 = IA2 * (seed2 - k * IQ2) - k * IR2;
    
    if (seed2 < 0) seed2 += IM2;
    
    j = long(iy / NDIV);
    iy = iv[j] - seed2;
    iv[j] = seed;
    
    if (iy < 1) iy += IMM1;
    
    if ((temp = double(AM * iy)) > RNMX)
    {
        return (double)RNMX;
    }
    else
    {
        return (double)temp;
    }
}
/* ------------------------------------------------------------------------------- */

double RandNum::pcg_uniform()
{
    return dist(rng);
}
/* ------------------------------------------------------------------------------- */

double RandNum::gaussian()
{
    double first,v1,v2,rsq,fac;
    
    if (!save) {
        
        int again = 1;
        
        while (again) {
            
            v1 = 2.0 * uniform() - 1.0;
            v2 = 2.0 * uniform() - 1.0;
            rsq = v1 * v1 + v2 * v2;
            
            if (rsq < 1.0 && rsq != 0.0) again = 0;
        }
        fac = std::sqrt(-2.0 * std::log(rsq) / rsq);
        second = v1 * fac;
        first = v2 * fac;
        save = 1;
        
    }
    else {
        first = second;
        save = 0;
    }
    
    return first;
}
/* ------------------------------------------------------------------------------- */

VectorXd RandNum::sphere()
{
    VectorXd n = VectorXd::Zero(3);
    
    double beta_1, beta_2;
    double s1 = 10.0;
    
    do
    {
        beta_1 = 2.0 * uniform() - 1.0;
        beta_2 = 2.0 * uniform() - 1.0;
        s1 = beta_1 * beta_1 + beta_2 * beta_2;
        
    } while (s1 > 1.0);
    
    double s2 = 2.0 * sqrt(1.0 - s1);
    
    n[0] = beta_1 * s2;
    n[1] = beta_2 * s2;
    n[2] = 1 - 2.0 * s1;
    
    return n;
}

/* ------------------------------------------------------------------------------- */

long RandNum::getRandState()
{
    if (parameters->rngType == RAN2)  {
        
        return rSeed;
    }
    else {
        return 0;
    }
}
/* ------------------------------------------------------------------------------- */


// ==============================================================================
//
//  Scatter2d.cpp
//  QTR
//
//  Created by Albert Lu on 4/21/26.
//  alu@tacc.utexas.edu
//
//  Last modified on 9/13/26
//
//  Note:
//
// ==============================================================================

#include <algorithm>
#include <cctype>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <new>
#include <omp.h>
#if defined(__GNUC__) && !defined(__clang__)
#include <parallel/algorithm>
#define QTR_SORT __gnu_parallel::sort
#else
#define QTR_SORT std::sort
#endif
#include <vector>
#include <string>

#include "Constants.h"
#include "Containers.h"
#include "Error.h"
#include "Log.h"
#include "Parameters.h"
#include "Scatter2d.h"

using namespace QTR_NS;
using std::vector;
using std::complex;
using std::real;
using std::imag;
using std::max;
using std::min;
using std::nothrow;

#define BIG_NUMBER 2147483647

// iostream extraction treats underflow (for example 1e-308) as a hard
// parse failure on macOS.  The legacy Vista files contain such values, so
// parse tokens with strtod and deliberately accept ERANGE underflow.
static bool read_data_row(std::istream &input, double &real_value,
                          double &imag_value, bool require_imaginary)
{
    std::string line;
    while (std::getline(input, line)) {
        const char *cursor = line.c_str();
        while (*cursor && std::isspace(static_cast<unsigned char>(*cursor)))
            ++cursor;
        if (!*cursor)
            continue;

        char *end = nullptr;
        real_value = std::strtod(cursor, &end);
        if (end == cursor)
            return false;
        cursor = end;

        while (*cursor && std::isspace(static_cast<unsigned char>(*cursor)))
            ++cursor;
        if (*cursor) {
            imag_value = std::strtod(cursor, &end);
            if (end == cursor)
                return false;
        } else {
            imag_value = 0.0;
        }

        if (require_imaginary && !*cursor)
            return false;
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------------------- */

Scatter2d::Scatter2d(class QTR *q)
{
    qtr = q;
    err = qtr->error;
    log = qtr->log;
    parameters = qtr->parameters;
    init();
} 
/* ------------------------------------------------------------------------------- */

Scatter2d::~Scatter2d()
{     
    return;
}
/* ------------------------------------------------------------------------------- */

void Scatter2d::init()
{
    log->log("\n\n[Scatter2d] INIT starts ...\n");
    log->log("\n\n[Scatter2d] Potential type: custom-defined in PES.dat\n");

    // General parameters
    I = {0,1}; // sqrt(-1)
    xZERO = {0,0}; // complex zero
    DIMENSIONS = parameters->scxd_dimensions;
    if (DIMENSIONS != 2) {
        log->log("[Scatter2d] ERROR: only dimensions=2 is supported by Scatter2d.\n");
        err->abort_all();
    }
    EDGE = parameters->scxd_edge > 2 ? parameters->scxd_edge : 2;
    PERIOD = parameters->scxd_period;
    PRINT_PERIOD = parameters->scxd_printperiod;
    TIME = parameters->scxd_Tf;
    kk = parameters->scxd_k;
    QUIET = parameters->quiet;
    TIMING = parameters->timing;
    DEBUG = parameters->scxd_isDebug;
    isToggle = parameters->scxd_isToggle;
    isReNorm = parameters->scxd_isReNorm;
    isFluxNProb = parameters->scxd_isFluxNProb;
    isFlux = parameters->scxd_isFlux;
    isCondProb = parameters->scxd_isCondProb;
    isPrintEdge = parameters->scxd_isPrintEdge;
    isPrintDensity = parameters->scxd_isPrintDensity;
    isPrintWaveFunc = parameters->scxd_isPrintWavefunc;

    if (PERIOD <= 0 || PRINT_PERIOD <= 0 || kk <= 0.0 || TIME < 0.0 ||
        parameters->scxd_h1 <= 0.0 || parameters->scxd_h2 <= 0.0 ||
        parameters->scxd_m1 <= 0.0 || parameters->scxd_m2 <= 0.0 ||
        parameters->scxd_hb <= 0.0) {
        log->log("[Scatter2d] ERROR: invalid time-step, period, grid-spacing, or mass.\n");
        err->abort_all();
    }

    log->log("[Scatter2d] DIMENSIONS: %d\n", DIMENSIONS);
    log->log("[Scatter2d] EDGE: %d\n", EDGE);
    log->log("[Scatter2d] PERIOD: %d\n", PERIOD);
    log->log("[Scatter2d] PRINT_PERIOD: %d\n", PRINT_PERIOD);
      
    // Grid size
    H.resize(DIMENSIONS);
    Hi.resize(DIMENSIONS);
    Hisq.resize(DIMENSIONS);
    H[0] = parameters->scxd_h1;
    H[1] = parameters->scxd_h2;

    for (int i = 0; i < DIMENSIONS; i ++)  {
        Hi[i] = 1.0 / H[i];
        Hisq[i] = 1.0 / (H[i] * H[i]);
    }

    // Domain size and # grids
    Box.resize(DIMENSIONS * 2);
    Box[0] = parameters->scxd_xi1;
    Box[1] = parameters->scxd_xf1;
    Box[2] = parameters->scxd_xi2;
    Box[3] = parameters->scxd_xf2;

    BoxShape.resize(DIMENSIONS);
    GRIDS_TOT = 1;
    log->log("[Scatter2d] Number of grids = (");

    for (int i = 0; i < DIMENSIONS; i ++)  {

        BoxShape[i] = (int)std::round((Box[2 * i + 1] - Box[2 * i]) / H[i]) + 1;
        GRIDS_TOT *= BoxShape[i];

        if ( i < DIMENSIONS - 1 )
            log->log("%d, ", BoxShape[i]);
        else
            log->log("%d)\n", BoxShape[i]);
    }

    if (DIMENSIONS==2)  {
        O1 = BoxShape[0] * BoxShape[1];
        M1 = BoxShape[1];
        W1 = BoxShape[1];
    }

    if (BoxShape[0] <= 2 * EDGE + 1 || BoxShape[1] <= 2 * EDGE + 1) {
        log->log("[Scatter2d] ERROR: grid is too small for the requested edge/stencil width.\n");
        err->abort_all();
    }

    // Parameters
    hb = parameters->scxd_hb;
    m1 = parameters->scxd_m1;
    m2 = parameters->scxd_m2;
    log->log("[Scatter2d] hb: %lf\n", hb);
    log->log("[Scatter2d] m1: %lf\n", m1);
    log->log("[Scatter2d] m2: %lf\n", m2);

    // Truncate parameters
    isFullGrid = parameters->scxd_isFullGrid;
    TolH = parameters->scxd_TolH;    // Tolerance of probability density for Zero point Cutoff
    TolL = parameters->scxd_TolL;    // Tolerance of probability density for Edge point
    TolHd = parameters->scxd_TolHd;  // Tolerance of probability first diff for Zero point Cutoff
    TolLd = parameters->scxd_TolLd;  // Tolerance of probability density for Edge point
    ExReduce = parameters->scxd_ExReduce; //Extrapolation reduce factor
    ExLimit = parameters->scxd_ExLimit;   //Extrapolation counts limit

    if (ExLimit < 0) {
        log->log("[Scatter2d] ERROR: ExLimit must be non-negative.\n");
        err->abort_all();
    }

    log->log("[Scatter2d] isFullGrid: %d\n", (int)isFullGrid);
    log->log("[Scatter2d] isToggle: %d\n", (int)isToggle);

    if (isToggle)  {
        toggle_threshold = parameters->scxd_toggle_threshold;
        log->log("[Scatter2d] Toggle threshold: %lf\n", toggle_threshold);
    }
    log->log("[Scatter2d] TolH: %e\n", TolH);
    log->log("[Scatter2d] TolL: %e\n", TolL);
    log->log("[Scatter2d] TolHd: %e\n", TolHd);
    log->log("[Scatter2d] TolLd: %e\n", TolLd);
    log->log("[Scatter2d] ExReduce: %lf\n", ExReduce);
    log->log("[Scatter2d] ExLimit: %d\n", ExLimit);

    // Original physical observables: one dissociation and three diffusion
    // saddle points, plus the six original conditional-probability regions.
    idx_abs_1 = parameters->scxd_idx_abs_1;
    idx_dsc_1 = parameters->scxd_idx_dsc_1;
    idx_dsc_2 = parameters->scxd_idx_dsc_2;
    idx_df1_1 = parameters->scxd_idx_df1_1;
    idx_df1_2 = parameters->scxd_idx_df1_2;
    idx_df2_1 = parameters->scxd_idx_df2_1;
    idx_df2_2 = parameters->scxd_idx_df2_2;
    idx_df3_1 = parameters->scxd_idx_df3_1;
    idx_df3_2 = parameters->scxd_idx_df3_2;

    if (isFlux || isCondProb) {
        const bool validOriginalObservables =
            idx_abs_1 >= 0 && idx_abs_1 <= BoxShape[0] &&
            idx_dsc_1 >= 0 && idx_dsc_1 < BoxShape[0] &&
            idx_df1_1 >= 0 && idx_df1_1 < BoxShape[0] &&
            idx_df2_1 >= 0 && idx_df2_1 < BoxShape[0] &&
            idx_df3_1 >= 0 && idx_df3_1 < BoxShape[0] &&
            idx_dsc_2 >= 2 && idx_dsc_2 + 2 < BoxShape[1] &&
            idx_df1_2 >= 2 && idx_df1_2 + 2 < BoxShape[1] &&
            idx_df2_2 >= 2 && idx_df2_2 + 2 < BoxShape[1] &&
            idx_df3_2 >= 2 && idx_df3_2 + 2 < BoxShape[1] &&
            idx_dsc_2 <= idx_df1_2 && idx_df1_2 <= idx_df2_2 &&
            idx_df2_2 <= idx_df3_2;
        if (!validOriginalObservables) {
            log->log("[Scatter2d] ERROR: invalid original flux/probability indices.\n");
            err->abort_all();
        }
    }

    log->log("[Scatter2d] idx_abs_1: %d\n", idx_abs_1);
    log->log("[Scatter2d] (idx_dsc_1,idx_dsc_2): %d %d\n", idx_dsc_1, idx_dsc_2);
    log->log("[Scatter2d] (idx_df1_1,idx_df1_2): %d %d\n", idx_df1_1, idx_df1_2);
    log->log("[Scatter2d] (idx_df2_1,idx_df2_2): %d %d\n", idx_df2_1, idx_df2_2);
    log->log("[Scatter2d] (idx_df3_1,idx_df3_2): %d %d\n", idx_df3_1, idx_df3_2);

    // Target-only desorption/dissociation planes remain available for the
    // optional re-normalization path below.
    idx_des_2_min = parameters->scxd_idx_des_2_min;
    idx_des_2_max = parameters->scxd_idx_des_2_max;
    idx_des_1 = parameters->scxd_idx_des_1;
    idx_dis_1_min = parameters->scxd_idx_dis_1_min;
    idx_dis_1_max = parameters->scxd_idx_dis_1_max;
    idx_dis_2 = parameters->scxd_idx_dis_2;

    if (isReNorm) {
        const bool validFluxPlane =
            idx_des_1 >= 2 && idx_des_1 + 2 < BoxShape[0] &&
            idx_dis_2 >= 2 && idx_dis_2 + 2 < BoxShape[1] &&
            idx_des_2_min >= 0 && idx_des_2_min < idx_des_2_max &&
            idx_des_2_max <= BoxShape[1] &&
            idx_dis_1_min >= 0 && idx_dis_1_min < idx_dis_1_max &&
            idx_dis_1_max <= BoxShape[0];
        if (!validFluxPlane) {
            log->log("[Scatter2d] ERROR: invalid flux-plane indices or ranges.\n");
            err->abort_all();
        }
    }
    log->log("[Scatter2d] idx_des_1: %d\n", idx_des_1);
    log->log("[Scatter2d] (idx_des_2_min,idx_des_2_max): %d %d\n", idx_des_2_min, idx_des_2_max);
    log->log("[Scatter2d] idx_dis_2: %d\n", idx_dis_2);
    log->log("[Scatter2d] (idx_dis_1_min,idx_dis_1_max): %d %d\n", idx_dis_1_min, idx_dis_1_max);

    log->log("[Scatter2d] INIT done.\n\n");
}
/* ------------------------------------------------------------------------------- */

void Scatter2d::Evolve()
{
    #pragma omp declare reduction (merge : MeshIndex : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

    log->log("[Scatter2d] Evolve starts ...\n");

    // Files
    FILE *pfile;

    // Variables 
    int count;
    int n1, n2;
    int ta_size, tb_size;
    int x1_max, x1_min, x2_max, x2_min;
    int x1_max_tmp, x1_min_tmp, x2_max_tmp, x2_min_tmp;
    int x1_flux_min, x1_flux_max;
    int x2_flux_min, x2_flux_max;
    double D1 = 1.0 / 12.0;
    double D2 = 2.0 / 3.0;
    double KC1 = 1.0 / 12.0;
    double KC2 = 4.0 / 3.0;
    double KC3 = 5.0 / 2.0;
    complex<double> sum;

    double kh2m1 = kk * hb / 2 / m1;
    double kh2m2 = kk * hb / 2 / m2;
    double kh2m1hisq0 = kh2m1 * Hisq[0];
    double kh2m2hisq1 = kh2m2 * Hisq[1];
    double k2hb = kk / hb;
    double h2m1h0 = hb / m1 * Hi[0];
    double h2m2h1 = hb / m2 * Hi[1];
    double flux_des;
    double flux_dis;
    double prob_remain = 1.0;
    double norm; 

    bool b1, b2, b3, b4, b5;
    bool isEmpty;
    
    // Timing variables
    double t_0_begin, t_1_begin, t_2_begin;    
    double t_0_end, t_1_end, t_2_end;
    double t_0_elapsed = 0.0;
    double t_1_elapsed = 0.0;
    double t_2_elapsed = 0.0;

    // Core computation time (RK4, normalization, initialization, etc)
    double t_full = 0.0;
    double t_truncate = 0.0;

    // Overhead time (truncate)
    double t_overhead = 0.0;

    // Constants
    double TolHd2H0, TolHd2H1;
    double TolLd2H0, TolLd2H1;

    // temporary index container
    MeshIndex tmpVec; 

    // Boundary layer container for extrapolation loop
    MeshIndex ExBD, ExBD_A;     
 
    //  1d Grid vector and indices
    VectorXi grid;
    int g1,g2;
    double xx1,xx2;
    std::complex<double> d2f1, d2f2, d1f1, d1f2;
    complex<double> f0;

    // Vector iterater
    vector<int>::iterator it;

    // Extrapolation 
    int min_dir;
    int Excount;
    bool isExtrapolate; 
    bool isFirstExtrp;
    complex<double> val, val_min;
    double val_min_abs;
    vector<complex<double>> ExTBL;

    // Neighborlist
    int nneigh = 0;
    vector<vector<int>> neighlist;
    vector<int> neighs(DIMENSIONS);
    log->log("[Scatter2d] Initializing containers ...\n");

    // Initialize containers

    t_0_begin = omp_get_wtime();

    bool *TAMask = nullptr;

    if ( !isFullGrid )  {
        TAMask = new(nothrow) bool[O1]();
    }

    double *PF = new(nothrow) double[O1];
    if (!PF) log->log("[Scatter2d] Fail allocating array PF\n");

    complex<double> *F = new(nothrow) complex<double>[O1];
    if (!F) log->log("[Scatter2d] Fail allocating array F\n");

    complex<double> *POT = new(nothrow) complex<double>[O1];
    if (!POT) log->log("[Scatter2d] Fail allocating array POT\n");

    complex<double> *FF = new(nothrow) complex<double>[O1];
    if (!FF) log->log("[Scatter2d] Fail allocating array FF\n");

    complex<double> *KK1 = new(nothrow) complex<double>[O1];
    if (!KK1) log->log("[Scatter2d] Fail allocating array KK1\n");

    complex<double> *KK2 = new(nothrow) complex<double>[O1];
    if (!KK2) log->log("[Scatter2d] Fail allocating array KK2\n");

    complex<double> *KK3 = new(nothrow) complex<double>[O1];
    if (!KK3) log->log("[Scatter2d] Fail allocating array KK3\n");

    complex<double> *KK4 = new(nothrow) complex<double>[O1];
    if (!KK4) log->log("[Scatter2d] Fail allocating array KK4\n");

    if ((!isFullGrid && !TAMask) || !PF || !F || !POT || !FF ||
        !KK1 || !KK2 || !KK3 || !KK4) {
        log->log("[Scatter2d] ERROR: unable to allocate simulation arrays.\n");
        err->abort_all();
    }

    #pragma omp parallel for schedule(static,4)
    for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
        for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {
            PF[i1*W1+i2] = 0.0;
            F[i1*W1+i2] = xZERO;
            POT[i1*W1+i2] = xZERO;
            FF[i1*W1+i2] = xZERO;
            KK1[i1*W1+i2] = xZERO;
            KK2[i1*W1+i2] = xZERO;
            KK3[i1*W1+i2] = xZERO;
            KK4[i1*W1+i2] = xZERO;
        }
    }

    if ( !isFullGrid )  {

        t_1_begin = omp_get_wtime();
        
        #pragma omp parallel for schedule(static,4)
        for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
            for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {
                TAMask[i1*W1+i2] = 0;
            }
        }

        nneigh = 0;

        for (int d = 1; d <= 4; d ++)  {
            for (n1 = -d; n1 <= d; n1 ++)  {

                n2 = d - abs(n1);

                if (n2 != 0)  {
                    neighs = {n1,n2};
                    neighlist.push_back(neighs);
                    neighs = {n1,-n2};
                    neighlist.push_back(neighs);
                    nneigh += 2;
                }
                else  {
                    neighs = {n1,0};
                    neighlist.push_back(neighs);
                    nneigh += 1;
                }
            }
        }
        log->log("[Scatter2d] nneigh = %d\n", nneigh); 
        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        t_overhead += t_1_elapsed;
    }
    t_0_end = omp_get_wtime();
    t_0_elapsed = t_0_end - t_0_begin;
    t_full += t_0_elapsed;

    if ( !isFullGrid )
        t_truncate += t_0_elapsed - t_1_elapsed; // subtract overhead

    if (!QUIET && TIMING) log->log("[Scatter2d] Elapsed time (initializing containers) = %.4e sec\n\n", t_0_elapsed); 

    // .........................................................................................

    log->log("[Scatter2d] Initializing custom-defined wavefunction from IniCndt.dat...\n");  

    t_1_begin = omp_get_wtime();

    // Initialize wavefunction

    std::ifstream wfn_val("IniCndt.dat");
    if (!wfn_val) {
        log->log("[Scatter2d] ERROR: cannot open IniCndt.dat.\n");
        err->abort_all();
    }
    
    double wfn_real, wfn_imgn;

    for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
        for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {

            if (!read_data_row(wfn_val, wfn_real, wfn_imgn, true)) {
                log->log("[Scatter2d] ERROR: IniCndt.dat read failed at grid (%d,%d), index %d of %d.\n", i1, i2, i1*W1+i2, O1);
                err->abort_all();
            }
            
            if (i1 < EDGE || i1 > BoxShape[0]-EDGE-1 || i2 < EDGE || i2 > BoxShape[1]-EDGE-1)
                F[i1*W1+i2] = xZERO;
            else
                F[i1*W1+i2] = {wfn_real, wfn_imgn};
        }
    }

    log->log("[Scatter2d] Initializing custom-defined potential energy surface from PES.dat...\n"); 

    // Initialize potential energy surface

    std::ifstream pot_val("PES.dat");
    if (!pot_val) {
        log->log("[Scatter2d] ERROR: cannot open PES.dat.\n");
        err->abort_all();
    }
    
    double pot_real, pot_imgn;
    
    for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
        for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {
    
            if (!read_data_row(pot_val, pot_real, pot_imgn, false)) {
                log->log("[Scatter2d] ERROR: invalid PES.dat row at grid index %d.\n", i1*W1+i2);
                err->abort_all();
            }
            POT[i1*W1+i2] = {pot_real, pot_imgn};
        }
    }

    // Normalization
    norm = 0.0;

    #pragma omp parallel for reduction (+:norm) schedule(static,4)
    for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
        for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
            norm += std::norm(F[i1*W1+i2]);
        }
    }
    norm *= H[0] * H[1];
    log->log("[Scatter2d] Normalization factor N = %.16e\n",norm);

    /* ------------------------------- */

    if (norm <= 0.0 || !std::isfinite(norm)) {
        log->log("[Scatter2d] ERROR: initial wavefunction has zero or invalid norm.\n");
        err->abort_all();
    }
    norm = 1.0 / sqrt(norm);
    log->log("[Scatter2d] Normalization factor 1/sqrt(N) = %.16e\n",norm);

    #pragma omp parallel for schedule(static,4)
    for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
        for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
            F[i1*W1+i2] = norm * F[i1*W1+i2];
            PF[i1*W1+i2] = std::norm(F[i1*W1+i2]);
        }
    }

    t_1_end = omp_get_wtime();
    t_1_elapsed = t_1_end - t_1_begin;
    t_full += t_1_elapsed;
    t_truncate += t_1_elapsed;
    if (!QUIET && TIMING) log->log("[Scatter2d] Elapsed time (initializing wavefunction) = %.4e sec\n\n", t_1_elapsed); 

    // .........................................................................................

    // Initial truncation & edge point check

    if ( !isFullGrid )
    {
        t_2_begin = omp_get_wtime();

        log->log("[Scatter2d] Initial truncation ...\n");

        TolHd2H0 = TolHd * 2.0 * H[0];
        TolHd2H1 = TolHd * 2.0 * H[1];
        TolLd2H0 = TolLd * 2.0 * H[0];
        TolLd2H1 = TolLd * 2.0 * H[1];

        x1_min = BIG_NUMBER;
        x2_min = BIG_NUMBER;
        x1_max = -BIG_NUMBER;
        x2_max = -BIG_NUMBER;

        // Truncation

        ta_size = 0;

        t_1_begin = omp_get_wtime();

        #pragma omp parallel for reduction(+: ta_size) schedule(static,4)
        for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
            for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                if (PF[i1*W1+i2] < TolH &&
                    std::abs(F[(i1+1)*W1+i2] - F[(i1-1)*W1+i2]) < TolHd2H0 &&
                    std::abs(F[i1*W1+(i2+1)] - F[i1*W1+(i2-1)]) < TolHd2H1) {
                    PF[i1*W1+i2] = 0.0;
                    continue;
                }
                TAMask[i1*W1+i2] = 1;
                ta_size += 1;
            }
        }

        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        if (!QUIET && TIMING) log->log("[Scatter2d] Elapsed time (initializing truncation A-1) = %.4e sec\n\n", t_1_elapsed); 

        // TA box and TB

        if (ta_size == 0)  {
            tb_size = 0;
            log->log("[Scatter2d] TA is empty\n");
        }
        else  {
            t_1_begin = omp_get_wtime();

            #pragma omp parallel for reduction(min: x1_min,x2_min) reduction(max: x1_max,x2_max) schedule(static,4)
            for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
                for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                    if (TAMask[i1*W1+i2])  {
                        if (i1 < x1_min)  x1_min = i1;
                        if (i1 > x1_max)  x1_max = i1;
                        if (i2 < x2_min)  x2_min = i2;
                        if (i2 > x2_max)  x2_max = i2;
                    }
                    else  {
                        F[i1*W1+i2] = xZERO;
                    }
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            if (!QUIET && TIMING) log->log("[Scatter2d] Elapsed time (initializing truncation A-2) = %.4e sec\n\n", t_1_elapsed); 
            t_1_begin = omp_get_wtime();

            // TB
            #pragma omp parallel for reduction(merge: tmpVec) schedule(static,4)
            for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                    if (TAMask[i1*W1+i2]) {
                        if( \
                            !TAMask[(i1+1)*W1+i2] || !TAMask[(i1-1)*W1+i2] || \
                            !TAMask[i1*W1+(i2+1)] || !TAMask[i1*W1+(i2-1)])
                            tmpVec.push_back(i1*W1+i2);
                    }
                }
            }
            tmpVec.swap(TB);
            tmpVec.clear();
            tb_size = TB.size();

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            if (!QUIET && TIMING) log->log("[Scatter2d] Elapsed time (initializing truncation A-3) = %.4e sec\n\n", t_1_elapsed); 
            t_1_begin = omp_get_wtime();

            // TA expansion

            #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2) schedule(static,4)
            for (int i = 0; i < TB.size(); i++)  {
                g1 = (int)(TB[i] / M1);
                g2 = (int)(TB[i] % M1);

                if (g1+1 < BoxShape[0]-EDGE-1 && !TAMask[(g1+1)*W1+g2])
                    tmpVec.push_back(GridToIdx(g1+1,g2));
                if (g1-1 > EDGE && !TAMask[(g1-1)*W1+g2])
                    tmpVec.push_back(GridToIdx(g1-1,g2));
                if (g2+1 < BoxShape[1]-EDGE-1 && !TAMask[g1*W1+(g2+1)])
                    tmpVec.push_back(GridToIdx(g1,g2+1));
                if (g2-1 > EDGE && !TAMask[g1*W1+(g2-1)])
                    tmpVec.push_back(GridToIdx(g1,g2-1));
            }

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            if (!QUIET && TIMING) log->log("[Scatter2d] Elapsed time (initializing truncation A-4) = %.4e sec\n\n", t_1_elapsed); 
            t_1_begin = omp_get_wtime();

            // Find unique elements
            QTR_SORT(tmpVec.begin(),tmpVec.end());
            it = std::unique (tmpVec.begin(), tmpVec.end()); 
            tmpVec.resize(std::distance(tmpVec.begin(),it));

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            if (!QUIET && TIMING) log->log("[Scatter2d] Elapsed time (initializing truncation A-5) = %.4e sec\n\n", t_1_elapsed);
            t_1_begin = omp_get_wtime();

            // Update TA box

            #pragma omp parallel for reduction(min: x1_min,x2_min) \
                                     reduction(max: x1_max,x2_max) \
                                     private(g1,g2) schedule(static,4)
            for (int i = 0; i < tmpVec.size(); i ++)  {
                g1 = (int)(tmpVec[i] / M1);
                g2 = (int)(tmpVec[i] % M1);

                if (!TAMask[g1*W1+g2]) {
                    TAMask[g1*W1+g2] = 1;
                    PF[g1*W1+g2] = 0.0;
                }

                x1_min = (g1 < x1_min) ? g1 : x1_min;
                x2_min = (g2 < x2_min) ? g2 : x2_min;
                x1_max = (g1 > x1_max) ? g1 : x1_max;
                x2_max = (g2 > x2_max) ? g2 : x2_max;

                if (DEBUG)
                    log->log("[DEBUG-TG3-A-%d] %.12e %.12e\n",0,Box[0]+g1*H[0],Box[2]+g2*H[1]);
            }
            tmpVec.clear();

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            if (!QUIET && TIMING) log->log("[Scatter2d] Elapsed time (initializing truncation A-6) = %.4e sec\n\n", t_1_elapsed); 
            t_1_begin = omp_get_wtime();

            // Update ta_size

            ta_size = 0;

            #pragma omp parallel for reduction(+: ta_size) schedule(static,4)
            for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                    if (TAMask[i1*W1+i2]) {
                        ta_size += 1;
                    }
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            if (!QUIET && TIMING) log->log("[Scatter2d] Elapsed time (initializing truncation A-7) = %.4e sec\n\n", t_1_elapsed); 
        }
        log->log("[Scatter2d] TA size = %d, TB size = %d\n", ta_size, tb_size);

        if (ta_size != 0)
            log->log("[Scatter2d] TA Range [%d, %d][%d, %d]\n", x1_min, x1_max, x2_min, x2_max);

        if (ta_size == 0)  {
            log->log("[Scatter2d] STOP: TA is empty\n");
            std::exit(EXIT_FAILURE);
        }

        log->log("[Scatter2d] Range [%d, %d][%d, %d]\n", x1_min, x1_max, x2_min, x2_max);

        t_2_end = omp_get_wtime();
        t_2_elapsed = t_2_end - t_2_begin;
        t_overhead += t_2_elapsed;

        if (!QUIET && TIMING)  {
            log->log("[Scatter2d] Elapsed time (initial truncation) = %.4e sec\n\n", t_1_elapsed);
            log->log("[Scatter2d] Initialization core computation time: %.4e sec\n", t_truncate);            
            log->log("[Scatter2d] Initialization overhead: %.4e sec\n", t_overhead); 
        }
    }
    else  // Full grid approach
    {
        log->log("[Scatter2d] Initialization core computation time: %.4e sec\n", t_full);   
    }

    // Original physical observables.  Keep one scalar per sampling point in
    // the original file names so existing analysis scripts remain usable.
    FILE *flux_d_file = nullptr;
    FILE *flux_1_file = nullptr;
    FILE *flux_2_file = nullptr;
    FILE *flux_3_file = nullptr;
    FILE *prob_r_a_file = nullptr;
    FILE *prob_a_d_file = nullptr;
    FILE *prob_d_1_file = nullptr;
    FILE *prob_1_2_file = nullptr;
    FILE *prob_2_3_file = nullptr;
    FILE *prob_3_p_file = nullptr;
    const bool write_observable_files = (parameters->me == 0);

    auto open_observable_file = [&](const char *filename) -> FILE * {
        FILE *file = fopen(filename, "w");
        if (file == nullptr) {
            log->log("[Scatter2d] ERROR: cannot open %s for writing.\n", filename);
            err->abort_all();
        }
        return file;
    };

    if (write_observable_files && isFlux) {
        flux_d_file = open_observable_file("Flux_d.dat");
        flux_1_file = open_observable_file("Flux_1.dat");
        flux_2_file = open_observable_file("Flux_2.dat");
        flux_3_file = open_observable_file("Flux_3.dat");
    }
    if (write_observable_files && isCondProb) {
        prob_r_a_file = open_observable_file("Prob_r_a.dat");
        prob_a_d_file = open_observable_file("Prob_a_d.dat");
        prob_d_1_file = open_observable_file("Prob_d_1.dat");
        prob_1_2_file = open_observable_file("Prob_1_2.dat");
        prob_2_3_file = open_observable_file("Prob_2_3.dat");
        prob_3_p_file = open_observable_file("Prob_3_p.dat");
    }

    auto write_original_observables = [&](int step) {
        if (isFlux) {
            auto flux_at = [&](int i1, int i2) {
                const complex<double> derivative =
                    D1 * (F[i1*W1+i2-2] - F[i1*W1+i2+2])
                    - D2 * (F[i1*W1+i2-1] - F[i1*W1+i2+1]);
                return (hb / m2 / H[1]) *
                    imag(conj(F[i1*W1+i2]) * derivative);
            };

            const double flux_d = flux_at(idx_dsc_1, idx_dsc_2);
            const double flux_1 = flux_at(idx_df1_1, idx_df1_2);
            const double flux_2 = flux_at(idx_df2_1, idx_df2_2);
            const double flux_3 = flux_at(idx_df3_1, idx_df3_2);

            log->log("[Scatter2d] Flux_d = %e\n", flux_d);
            log->log("[Scatter2d] Flux_1 = %e\n", flux_1);
            log->log("[Scatter2d] Flux_2 = %e\n", flux_2);
            log->log("[Scatter2d] Flux_3 = %e\n", flux_3);

            if (write_observable_files) {
                fprintf(flux_d_file, "%.16e\n", flux_d);
                fprintf(flux_1_file, "%.16e\n", flux_1);
                fprintf(flux_2_file, "%.16e\n", flux_2);
                fprintf(flux_3_file, "%.16e\n", flux_3);
            }
        }

        if (isCondProb) {
            double prob_r_a = 0.0;
            double prob_a_d = 0.0;
            double prob_d_1 = 0.0;
            double prob_1_2 = 0.0;
            double prob_2_3 = 0.0;
            double prob_3_p = 0.0;

            #pragma omp parallel for reduction(+:prob_r_a) schedule(static,4)
            for (int i1 = idx_abs_1; i1 < BoxShape[0]; ++i1)
                for (int i2 = 0; i2 < BoxShape[1]; ++i2)
                    prob_r_a += PF[i1*W1+i2];

            #pragma omp parallel for reduction(+:prob_a_d) schedule(static,4)
            for (int i1 = 0; i1 < idx_abs_1; ++i1)
                for (int i2 = 0; i2 < idx_dsc_2; ++i2)
                    prob_a_d += PF[i1*W1+i2];

            #pragma omp parallel for reduction(+:prob_d_1) schedule(static,4)
            for (int i1 = 0; i1 < idx_abs_1; ++i1)
                for (int i2 = idx_dsc_2; i2 < idx_df1_2; ++i2)
                    prob_d_1 += PF[i1*W1+i2];

            #pragma omp parallel for reduction(+:prob_1_2) schedule(static,4)
            for (int i1 = 0; i1 < idx_abs_1; ++i1)
                for (int i2 = idx_df1_2; i2 < idx_df2_2; ++i2)
                    prob_1_2 += PF[i1*W1+i2];

            #pragma omp parallel for reduction(+:prob_2_3) schedule(static,4)
            for (int i1 = 0; i1 < idx_abs_1; ++i1)
                for (int i2 = idx_df2_2; i2 < idx_df3_2; ++i2)
                    prob_2_3 += PF[i1*W1+i2];

            #pragma omp parallel for reduction(+:prob_3_p) schedule(static,4)
            for (int i1 = 0; i1 < idx_abs_1; ++i1)
                for (int i2 = idx_df3_2; i2 < BoxShape[1]; ++i2)
                    prob_3_p += PF[i1*W1+i2];

            const double cell_area = H[0] * H[1];
            prob_r_a *= cell_area;
            prob_a_d *= cell_area;
            prob_d_1 *= cell_area;
            prob_1_2 *= cell_area;
            prob_2_3 *= cell_area;
            prob_3_p *= cell_area;

            log->log("[Scatter2d] Prob_r_a = %e\n", prob_r_a);
            log->log("[Scatter2d] Prob_a_d = %e\n", prob_a_d);
            log->log("[Scatter2d] Prob_d_1 = %e\n", prob_d_1);
            log->log("[Scatter2d] Prob_1_2 = %e\n", prob_1_2);
            log->log("[Scatter2d] Prob_2_3 = %e\n", prob_2_3);
            log->log("[Scatter2d] Prob_3_p = %e\n", prob_3_p);

            if (write_observable_files) {
                fprintf(prob_r_a_file, "%.16e\n", prob_r_a);
                fprintf(prob_a_d_file, "%.16e\n", prob_a_d);
                fprintf(prob_d_1_file, "%.16e\n", prob_d_1);
                fprintf(prob_1_2_file, "%.16e\n", prob_1_2);
                fprintf(prob_2_3_file, "%.16e\n", prob_2_3);
                fprintf(prob_3_p_file, "%.16e\n", prob_3_p);
            }
        }

        (void)step;
    };

    write_original_observables(0);
    // .........................................................................................

    // Time iteration 

    log->log("=======================================================\n\n"); 
    log->log("[Scatter2d] Time interation starts ...\n"); 
    log->log("[Scatter2d] Number of steps = %d\n\n", (int)(TIME / kk)); 
    log->log("=======================================================\n\n"); 

    for (int tt = 0; tt < (int)(TIME / kk); tt ++)
    {
        t_0_begin = omp_get_wtime(); 

        Excount = 0;

        if ( tt % PRINT_PERIOD == 0 )
        {
            if ( isPrintEdge  && !isFullGrid )  {

                pfile = fopen ("edge.dat","a");
                fprintf(pfile, "%d %lf %lu\n", tt, tt * kk, TB.size());

                for (int i = 0; i < TB.size(); i++)
                {
                    g1 = (int)(TB[i] / M1);
                    g2 = (int)(TB[i] % M1);
                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];
                    fprintf(pfile, "%d %d %lf %lf\n", g1,g2,xx1,xx2);            
                }
                fclose(pfile);
            }
            if ( isPrintDensity && !isFullGrid )  {

                pfile = fopen ("density.dat","a");
                fprintf(pfile, "%d %lf %d\n", tt, tt * kk, ta_size);

                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        if (TAMask[i1*W1+i2])  {
                            xx1 = Box[0] + i1 * H[0];
                            xx2 = Box[2] + i2 * H[1];
                            fprintf(pfile, "%d %d %lf %lf %.16e\n", i1, i2, xx1, xx2, PF[i1*W1+i2]);
                        }
                    }
                }
                fclose(pfile);
            }
            if ( isPrintDensity && isFullGrid )  {

                pfile = fopen ("density.dat","a");
                fprintf(pfile, "%d %lf %d\n", tt, tt * kk, O1);

                for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
                    for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {
                        xx1 = Box[0] + i1 * H[0];
                        xx2 = Box[2] + i2 * H[1];
                        fprintf(pfile, "%d %d %lf %lf %.16e\n", i1, i2, xx1, xx2, PF[i1*W1+i2]);
                    }
                }
                fclose(pfile);
            }

            if ( isPrintWaveFunc && !isFullGrid )  {

                pfile = fopen ("wave.dat","a");
                fprintf(pfile, "%d %lf %d\n", tt, tt * kk, ta_size);

                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        if (TAMask[i1*W1+i2])  {
                            xx1 = Box[0] + i1 * H[0];
                            xx2 = Box[2] + i2 * H[1];
                            fprintf(pfile, "%d %d %lf %lf %.16e %.16e\n", i1, i2, xx1, xx2, real(F[i1*W1+i2]), imag(F[i1*W1+i2]));
                        }
                    }
                }
                fclose(pfile);
            }
            if ( isPrintWaveFunc && isFullGrid && 0)  {

                pfile = fopen ("wave.dat","a");

                fprintf(pfile, "%d %lf %d\n", tt, tt * kk, O1);

                for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
                    for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {
                            xx1 = Box[0] + i1 * H[0];
                            xx2 = Box[2] + i2 * H[1];
                            fprintf(pfile, "%d %d %lf %lf %.16e %.16e\n", i1, i2, xx1, xx2, real(F[i1*W1+i2]), imag(F[i1*W1+i2]));
                    }
                }
                fclose(pfile);
            }
        }

        //  Computation of Flux and Probability

        // Check if TB of f is higher than TolL
        
        if ( !isFullGrid )
        {
            t_1_begin = omp_get_wtime();
            t_truncate = 0.0;
            t_overhead = 0.0;

            // TBL = Index of Extrapolating Edge points.
            // TBL_P = Index history of TBL of this iteration. To Prevent from extrapolating the same points multiple times.

            TBL.clear();

            const int nthreads = omp_get_max_threads();

            std::vector<std::vector<int>> localVecs(nthreads);

            #pragma omp parallel 
            {
                const int tid = omp_get_thread_num();

                auto &local = localVecs[tid];
                local.reserve(TB.size() / nthreads);

                #pragma omp for private(g1,g2,b1,b2,b3,b4,b5) schedule(static,4)
                for (int i = 0; i < TB.size(); i++)
                {
                    const int idx = TB[i];
                    g1 = idx / M1;
                    g2 = idx - g1 * M1;        

                    int center = g1*W1 + g2;

                    b1 = PF[center] >= TolL;
                    b2 = std::abs(F[center+W1] - F[center-W1]) >= TolLd2H0;
                    b3 = std::abs(F[center+1] - F[center-1]) >= TolLd2H1;
                    b4 = g1 > EDGE && g2 > EDGE;
                    b5 = (g1 < BoxShape[0]-EDGE-1) && (g2 < BoxShape[1]-EDGE-1);

                    bool keep = (b1||b2||b3) && b4 && b5;

                    if (keep)
                        local.push_back(TB[i]);
                }
            }

            size_t totalSize = 0;

            for (const auto &v : localVecs)
            {
                totalSize += v.size();
            }

            TBL.reserve(totalSize);

            // merge all thread-local vectors
            for (auto &v : localVecs)
            {
                TBL.insert(TBL.end(), v.begin(), v.end());
            }

            TBL_P = TBL;           

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-a-1: TBL) = %.4e sec\n", t_1_elapsed);   
        }
        else  
        {
            t_full = 0.0;
        }
        isExtrapolate = false;
        isFirstExtrp = true;
 
        // CASE 1: Truncating with extrapolation

        while ( !isFullGrid && TBL.size() != 0 && Excount < ExLimit )
        {
            // Extrapolation

            if ( TBL.size() != 0 )  {

                t_1_begin = omp_get_wtime();

                isExtrapolate = true;

                // Avoid unexpected arrangement of TBL
                QTR_SORT(TBL.begin(),TBL.end());
                it = std::unique (TBL.begin(), TBL.end()); 
                TBL.resize(std::distance(TBL.begin(),it)); 

                // Find extrapolation target
                // ExFF: Index of Extrapolated points
                ExFF.clear();

                //#pragma omp parallel for reduction(merge: tmpVec) private(g1,g2) schedule(static,4)
                for (int i = 0; i < TBL.size(); i++)  {

                    g1 = (int)(TBL[i] / M1);
                    g2 = (int)(TBL[i] % M1);

                    if ( g1-1 > EDGE && PF[(g1-1)*W1+g2] == 0.0 ) {
                        tmpVec.push_back(GridToIdx(g1-1,g2));
                    }
                    if ( g1+1 < BoxShape[0]-EDGE-1 && PF[(g1+1)*W1+g2] == 0.0 ) {
                        tmpVec.push_back(GridToIdx(g1+1,g2));
                    }
                    if ( g2-1 > EDGE && PF[g1*W1+(g2-1)] == 0.0 ) {
                        tmpVec.push_back(GridToIdx(g1,g2-1));
                    }
                    if ( g2+1 < BoxShape[1]-EDGE-1 && PF[g1*W1+(g2+1)] == 0.0 ) {
                        tmpVec.push_back(GridToIdx(g1,g2+1));
                    }
                }
                tmpVec.swap(ExFF);
                tmpVec.clear();

                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-b-1: ExFF A) = %.4e sec\n", t_1_elapsed);  

                if ( ExFF.size() > 0 )  {

                    t_1_begin = omp_get_wtime();

                    // ExFF & TBL set difference
                    tmpVec.resize(ExFF.size() + TBL.size());
                    QTR_SORT(TBL.begin(), TBL.end());
                    QTR_SORT(ExFF.begin(), ExFF.end());
                    it = std::set_difference(ExFF.begin(), ExFF.end(), TBL.begin(), TBL.end(), tmpVec.begin());
                    tmpVec.resize(it - tmpVec.begin()); 
                    tmpVec.swap(ExFF);
                    tmpVec.clear();

                    // Find unique elements
                    QTR_SORT(ExFF.begin(),ExFF.end());
                    it = std::unique (ExFF.begin(), ExFF.end()); 
                    ExFF.resize(std::distance(ExFF.begin(),it));

                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-b-2: ExFF A) = %.4e sec\n", t_1_elapsed);   

                    // Find the direction of Outer to Edge points
                    t_1_begin = omp_get_wtime();

                    Check.clear();
                    ExTBL.clear();

                    for (int i = 0; i < ExFF.size(); i ++)  {
                        Check.push_back(1);
                        ExTBL.push_back(xZERO);
                    }
                    
                    #pragma omp parallel for private(g1,g2,sum,count,isEmpty,val,val_min_abs,val_min,min_dir) schedule(static,4)
                    for (int i = 0; i < ExFF.size(); i ++)
                    {
                        g1 = (int)(ExFF[i] / M1);
                        g2 = (int)(ExFF[i] % M1);
                        sum = xZERO;
                        count = 0;
                        isEmpty = true;
                        val_min_abs = BIG_NUMBER;
                        val_min = BIG_NUMBER;
                        min_dir = -1;

                        if ( F[(g1-1)*W1+g2] != xZERO )  {

                            if ( std::abs(F[(g1-1)*W1+g2]) < val_min_abs &&  F[(g1-2)*W1+g2] != xZERO )  {
                                val_min_abs = std::abs(F[(g1-1)*W1+g2]);
                                val_min = F[(g1-1)*W1+g2];
                                min_dir = 0;
                            }
                            if ( F[(g1-2)*W1+g2] != xZERO )  {

                                val = exp( 2.0 * std::log(F[(g1-1)*W1+g2]) - std::log(F[(g1-2)*W1+g2]) );

                                if ( !(std::isnan(real(val)) || std::isnan(-real(val))) && !(std::isinf(real(val)) || std::isinf(-real(val))) && !(std::isnan(imag(val)) || std::isnan(-imag(val))) && !(std::isinf(imag(val)) || std::isinf(-imag(val))) )  
                                {
                                    sum += val;
                                    count += 1;
                                    isEmpty = false;
                                }
                            }
                        }

                        if ( F[(g1+1)*W1+g2] != xZERO )  {

                            if ( std::abs(F[(g1+1)*W1+g2]) < val_min_abs && F[(g1+2)*W1+g2] != xZERO )  {
                                val_min_abs = std::abs(F[(g1+1)*W1+g2]);
                                val_min = F[(g1+1)*W1+g2];
                                min_dir = 0;
                            }
                            if ( F[(g1+2)*W1+g2] != xZERO )  {

                                val = exp( 2.0 * std::log(F[(g1+1)*W1+g2]) - std::log(F[(g1+2)*W1+g2]) );

                                if ( !(std::isnan(real(val)) || std::isnan(-real(val))) && !(std::isinf(real(val)) || std::isinf(-real(val))) && !(std::isnan(imag(val)) || std::isnan(-imag(val))) && !(std::isinf(imag(val)) || std::isinf(-imag(val))) )  
                                {
                                    sum += val;
                                    count += 1;
                                    isEmpty = false;
                                }
                            }
                        }

                        if ( F[g1*W1+(g2-1)] != xZERO )  {

                            if ( std::abs(F[g1*W1+(g2-1)]) < val_min_abs && F[g1*W1+(g2-2)] != xZERO )  {
                                val_min_abs = std::abs(F[g1*W1+(g2-1)]);
                                val_min = F[g1*W1+(g2-1)];
                                min_dir = 1;
                            }
                            if ( F[g1*W1+(g2-2)] != xZERO )  {

                                val = exp( 2.0 * std::log(F[g1*W1+(g2-1)]) - std::log(F[g1*W1+(g2-2)]) );

                                if ( !(std::isnan(real(val)) || std::isnan(-real(val))) && !(std::isinf(real(val)) || std::isinf(-real(val))) && !(std::isnan(imag(val)) || std::isnan(-imag(val))) && !(std::isinf(imag(val)) || std::isinf(-imag(val))) )  
                                {
                                    sum += val;
                                    count += 1;
                                    isEmpty = false;
                                }
                            }
                        }

                        if ( F[g1*W1+(g2+1)] != xZERO )  {

                            if ( std::abs(F[g1*W1+(g2+1)]) < val_min_abs && F[g1*W1+(g2+2)] != xZERO )  {
                                val_min_abs = std::abs(F[g1*W1+(g2+1)]);
                                val_min = F[g1*W1+(g2+1)];
                                min_dir = 1;
                            }
                            if ( F[g1*W1+(g2+2)] != xZERO )  {

                                val = exp( 2.0 * std::log(F[g1*W1+(g2+1)]) - std::log(F[g1*W1+(g2+2)]) );

                                if ( !(std::isnan(real(val)) || std::isnan(-real(val))) && !(std::isinf(real(val)) || std::isinf(-real(val))) && !(std::isnan(imag(val)) || std::isnan(-imag(val))) && !(std::isinf(imag(val)) || std::isinf(-imag(val))) )  
                                {
                                    sum += val;
                                    count += 1;
                                    isEmpty = false;
                                }
                            }
                        }

                        if ( isEmpty )
                        {
                            Check[i] = 0; 

                        }  else  {

                            // Assume the probability of outer points always smaller than edge point.
                            // if larger, then set the edge point with smallest P to outer
                            // point instead of using the extrapolation result.
                            if ( std::abs(sum)/count > val_min_abs )  {
                                ExTBL[i] = val_min * exp(-ExReduce * H[min_dir]);
                            }
                            else  {
                                ExTBL[i] = sum / (double)count;
                            }
                        }
                    }

                    count = 0;

                    #pragma omp parallel for reduction (+:count) schedule(static,4)
                    for ( int i = 0; i < ExFF.size(); i++ )  {
                        if (Check[i] == 1)  {
                            F[ExFF[i]] = ExTBL[i];
                            count += 1;
                        }
                    }

                    if (count == 0)  {
                        ExFF.clear();
                        ExTBL.clear();
                        Check.clear();
                    }
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-b-3: ExFF A) = %.4e sec\n", t_1_elapsed);  
                }  // ExFF.size() > 0
            }  // if TBL.size() != 0

            // ............................................................................................. Extrapolation

            if ( isFirstExtrp )  {

                // Check Extending nonzero Area

                if ( ExFF.size() > 0 )  {

                    t_1_begin = omp_get_wtime();

                    tmpVec.clear();

                    #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2) schedule(static,4)
                    for (int i = 0; i < ExFF.size(); i++)  {
                        if (Check[i] == 1)  
                        {
                            g1 = (int)(ExFF[i] / M1);
                            g2 = (int)(ExFF[i] % M1);

                            if ( !TAMask[g1*W1+g2] )
                                tmpVec.push_back(g1*W1+g2);
                            if ( g1+1 < BoxShape[0]-EDGE-1 && !TAMask[(g1+1)*W1+g2] )
                                tmpVec.push_back((g1+1)*W1+g2);
                            if ( g1-1 > EDGE && !TAMask[(g1-1)*W1+g2] )
                                tmpVec.push_back((g1-1)*W1+g2);
                            if ( g2+1 < BoxShape[1]-EDGE-1 && !TAMask[g1*W1+(g2+1)] )
                                tmpVec.push_back(g1*W1+(g2+1));
                            if ( g2-1 > EDGE && !TAMask[g1*W1+(g2-1)] )
                                tmpVec.push_back(g1*W1+(g2-1));
                        }
                    }

                    for (int i = 0; i < tmpVec.size(); i ++)  {
                        g1 = (int)(tmpVec[i] / M1);
                        g2 = (int)(tmpVec[i] % M1);
                        TAMask[tmpVec[i]] = 1;
                        PF[tmpVec[i]] = 0.0;
                        F[tmpVec[i]] = xZERO;
                    
                        x1_min = (g1 < x1_min) ? g1 : x1_min;
                        x1_max = (g1 > x1_max) ? g1 : x1_max;
                        x2_min = (g2 < x2_min) ? g2 : x2_min;
                        x2_max = (g2 > x2_max) ? g2 : x2_max;
                    }
                    tmpVec.clear();

                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-c-1: CASE 1 TA A) = %.4e sec\n", t_1_elapsed);
                }

                if (DEBUG)  {
                    log->log("[DEBUG-TA-1A] [%d, %d][%d, %d]\n", x1_min, x1_max, x2_min, x2_max);
                }

                // Runge–Kutta 4
                #pragma omp parallel 
                {
                    #pragma omp single nowait
                    {
                        t_1_begin = omp_get_wtime();
                    }

                    // RK4-1
                    #pragma omp for private(d2f1,d2f2) schedule(static,4)
                    for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                        for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                            if (TAMask[i1*W1+i2])  {
                                d2f1 = (i1 >= EDGE && i1 < BoxShape[0] - EDGE) ? 
                                - KC1 * (F[(i1-2)*W1+i2] + F[(i1+2)*W1+i2])
                                + KC2 * (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2])
                                - KC3 * F[i1*W1+i2] : 
                                + F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] - 2.0 * F[i1*W1+i2];

                                d2f2 = (i2 >= EDGE && i2 < BoxShape[1] - EDGE) ? 
                                - KC1 * (F[i1*W1+(i2-2)] + F[i1*W1+(i2+2)])
                                + KC2 * (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)]) 
                                - KC3 * F[i1*W1+i2] : 
                                + F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] - 2.0 * F[i1*W1+i2];

                                KK1[i1*W1+i2] = I * (kh2m1hisq0 * d2f1 + kh2m2hisq1 * d2f2)  \
                                                - I * k2hb * POT[i1*W1+i2] * F[i1*W1+i2];

                                FF[i1*W1+i2] = F[i1*W1+i2] + KK1[i1*W1+i2] / 6.0;
                            }
                        }
                    }

                    #pragma omp single nowait
                    {
                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        t_truncate += t_1_elapsed;
                        if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-11: CASE 1 KK1) = %.4e sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();
                    }

                    // RK4-2
                    #pragma omp for private(d2f1, d2f2) schedule(static,4)
                    for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                        for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                            if (TAMask[i1*W1+i2])  {
                                d2f1 = (i1 >= EDGE && i1 < BoxShape[0] - EDGE) ?
                                    - KC1 * (F[(i1-2)*W1+i2] + F[(i1+2)*W1+i2] + 0.5 * (KK1[(i1-2)*W1+i2] + KK1[(i1+2)*W1+i2]))
                                    + KC2 * (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] + 0.5 * (KK1[(i1-1)*W1+i2] + KK1[(i1+1)*W1+i2]))
                                    - KC3 * (F[i1*W1+i2] + 0.5 * KK1[i1*W1+i2]) :
                                    + (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] + 0.5 * (KK1[(i1-1)*W1+i2] + KK1[(i1+1)*W1+i2])) 
                                    - 2.0 * F[i1*W1+i2] - KK1[i1*W1+i2];
                                
                                d2f2 = (i2 >= EDGE && i2 < BoxShape[1] - EDGE) ?
                                    - KC1 * (F[i1*W1+(i2-2)] + F[i1*W1+(i2+2)] + 0.5 * (KK1[i1*W1+(i2-2)] + KK1[i1*W1+(i2+2)]))
                                    + KC2 * (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] + 0.5 * (KK1[i1*W1+(i2-1)] + KK1[i1*W1+(i2+1)]))
                                    - KC3 * (F[i1*W1+i2] + 0.5 * KK1[i1*W1+i2]) :
                                    + (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] + 0.5 * (KK1[i1*W1+(i2-1)] + KK1[i1*W1+(i2+1)]))
                                    - 2.0 * F[i1*W1+i2] - KK1[i1*W1+i2];

                                KK2[i1*W1+i2] = I * (kh2m1hisq0 * d2f1 + kh2m2hisq1 * d2f2) \
                                                - I * k2hb * POT[i1*W1+i2] * (F[i1*W1+i2] + 0.5 * KK1[i1*W1+i2]);

                                FF[i1*W1+i2] += KK2[i1*W1+i2] / 3.0;
                            }
                        }
                    }

                    #pragma omp single nowait
                    {
                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        t_truncate += t_1_elapsed;
                        if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-12: CASE 1 KK2) = %.4e sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();
                    }

                    // RK4-3
                    #pragma omp for private(d2f1, d2f2) schedule(static,4)
                    for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                        for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                            if (TAMask[i1*W1+i2])  {
                                d2f1 = (i1 >= EDGE && i1<BoxShape[0] - EDGE) ?
                                    - KC1 * (F[(i1-2)*W1+i2] + F[(i1+2)*W1+i2] + 0.5 * (KK2[(i1-2)*W1+i2] + KK2[(i1+2)*W1+i2]))
                                    + KC2 * (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] + 0.5 * (KK2[(i1-1)*W1+i2] + KK2[(i1+1)*W1+i2]))
                                    - KC3 * (F[i1*W1+i2] + 0.5 * KK2[i1*W1+i2]) :
                                    + (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] + 0.5 * (KK2[(i1-1)*W1+i2] + KK2[(i1+1)*W1+i2]))
                                    - 2.0 * F[i1*W1+i2] - KK2[i1*W1+i2];

                                d2f2 = (i2 >= EDGE && i2 < BoxShape[1] - EDGE) ?
                                    - KC1 * (F[i1*W1+(i2-2)] + F[i1*W1+(i2+2)] + 0.5 * (KK2[i1*W1+(i2-2)] + KK2[i1*W1+(i2+2)]))
                                    + KC2 * (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] + 0.5 * (KK2[i1*W1+(i2-1)] + KK2[i1*W1+(i2+1)]))
                                    - KC3 * (F[i1*W1+i2] + 0.5 * KK2[i1*W1+i2]) :
                                    + (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] + 0.5 * (KK2[i1*W1+(i2-1)] + KK2[i1*W1+(i2+1)]))
                                    - 2.0 * F[i1*W1+i2] - KK2[i1*W1+i2];

                                KK3[i1*W1+i2] = I * (kh2m1hisq0 * d2f1 + kh2m2hisq1 * d2f2) \
                                                - I * k2hb * POT[i1*W1+i2] * (F[i1*W1+i2] + 0.5 * KK2[i1*W1+i2]);

                                FF[i1*W1+i2] += KK3[i1*W1+i2] / 3.0;
                            }
                        }
                    }

                    #pragma omp single nowait
                    {
                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        t_truncate += t_1_elapsed;
                        if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-13: CASE 1 KK3) = %.4e sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();
                    }

                    // RK4-4
                    #pragma omp for private(d2f1, d2f2) schedule(static,4)
                    for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                        for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                            if (TAMask[i1*W1+i2])  {
                                d2f1 = (i1 >= EDGE && i1 < BoxShape[0] - EDGE) ?
                                    - KC1 * (F[(i1-2)*W1+i2] + F[(i1+2)*W1+i2] + KK3[(i1-2)*W1+i2] + KK3[(i1+2)*W1+i2])
                                    + KC2 * (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] + KK3[(i1-1)*W1+i2] + KK3[(i1+1)*W1+i2])
                                    - KC3 * (F[i1*W1+i2] + KK3[i1*W1+i2]) :
                                    + (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] + KK3[(i1-1)*W1+i2] + KK3[(i1+1)*W1+i2])
                                    - 2.0 * (F[i1*W1+i2] + KK3[i1*W1+i2]);

                                d2f2 = (i2 >= EDGE && i2 < BoxShape[1] - EDGE) ?
                                    - KC1 * (F[i1*W1+(i2-2)] + F[i1*W1+(i2+2)] + KK3[i1*W1+(i2-2)] + KK3[i1*W1+(i2+2)])
                                    + KC2 * (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] + KK3[i1*W1+(i2-1)] + KK3[i1*W1+(i2+1)])
                                    - KC3 * (F[i1*W1+i2] + KK3[i1*W1+i2]) :
                                    + (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] + KK3[i1*W1+(i2-1)] + KK3[i1*W1+(i2+1)])
                                    - 2.0 * (F[i1*W1+i2] + KK3[i1*W1+i2]);

                                KK4[i1*W1+i2] = I * (kh2m1hisq0 * d2f1 + kh2m2hisq1 * d2f2) \
                                                - I * k2hb * POT[i1*W1+i2] * (F[i1*W1+i2] + KK3[i1*W1+i2]);

                                FF[i1*W1+i2] += KK4[i1*W1+i2] / 6.0;
                            }
                        }
                    }

                    #pragma omp single nowait
                    {
                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        t_truncate += t_1_elapsed;
                        if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-14: CASE 1 KK4) = %.4e sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();
                    }
                } // omp parallel

                isFirstExtrp = false;

            } // if ( isFirstExtrp)
            else if (ExFF.size() == 0)  {

                // In the case no valid ExFF found, reset TBL to break the while loop 
                TBL.clear();
            }
            else
            {
                // Extrapolation loop when multiple expanding occured

                if ( ExFF.size() > 0 )  {

                    t_1_begin = omp_get_wtime();

                    tmpVec.clear();

                    #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2) schedule(static,4)
                    for (int i = 0; i < ExFF.size(); i++)  {
                        if (Check[i] == 1)  
                        {
                            g1 = (int)(ExFF[i] / M1);
                            g2 = (int)(ExFF[i] % M1);

                            if ( !TAMask[g1*W1+g2] )
                                tmpVec.push_back(g1*W1+g2);
                            if ( g1+1 < BoxShape[0]-EDGE-1 && !TAMask[(g1+1)*W1+g2] )
                                tmpVec.push_back((g1+1)*W1+g2);
                            if ( g1-1 > EDGE && !TAMask[(g1-1)*W1+g2] )
                                tmpVec.push_back((g1-1)*W1+g2);
                            if ( g2+1 < BoxShape[1]-EDGE-1 && !TAMask[g1*W1+(g2+1)] )
                                tmpVec.push_back(g1*W1+(g2+1));
                            if ( g2-1 > EDGE && !TAMask[g1*W1+(g2-1)] )
                                tmpVec.push_back(g1*W1+(g2-1));
                        }
                    }

                    for (int i = 0; i < tmpVec.size(); i ++)  {

                        g1 = (int)(tmpVec[i] / M1); 
                        g2 = (int)(tmpVec[i] % M1);
                        TAMask[tmpVec[i]] = 1;

                        x1_min = (g1 < x1_min) ? g1 : x1_min;
                        x1_max = (g1 > x1_max) ? g1 : x1_max;
                        x2_min = (g2 < x2_min) ? g2 : x2_min;
                        x2_max = (g2 > x2_max) ? g2 : x2_max;
                    }
                    tmpVec.clear();
                    ExBD_A.clear();

                    #pragma omp parallel for reduction(merge: ExBD_A) private(g1,g2,n1,n2) schedule(static,4)
                    for (int i = 0; i < ExFF.size(); i++)
                    {
                        if (Check[i] == 1)  
                        {
                            g1 = (int)(ExFF[i] / M1);
                            g2 = (int)(ExFF[i] % M1);

                            ExBD_A.push_back(ExFF[i]);

                            for (int j = 0; j < nneigh; j ++)  {

                                n1 = neighlist[j][0];
                                n2 = neighlist[j][1];

                                if (TAMask[(g1+n1)*W1+(g2+n2)])
                                    ExBD_A.push_back((g1+n1)*W1+(g2+n2));
                            }
                        }
                    }

                    // Find unique elements (ExBD)
                    QTR_SORT(ExBD_A.begin(),ExBD_A.end());
                    it = std::unique (ExBD_A.begin(), ExBD_A.end()); 
                    ExBD_A.resize(std::distance(ExBD_A.begin(),it));

                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-cx-1: CASE 1 ExBD A) = %.4e sec\n", t_1_elapsed); 

                }  // if ExFF.size() > 0


                if ( ExFF.size() > 0 )  {

                    ExBD_A.swap(ExBD);

                    // Runge–Kutta 4

                    #pragma omp parallel
                    {
                        #pragma omp single nowait
                        {
                            t_1_begin = omp_get_wtime();
                        }
                        // RK4-1
                        #pragma omp for private(g1, g2, d2f1, d2f2) schedule(static,4)
                        for (int i = 0; i < ExBD.size(); i++)  {

                            g1 = (int)(ExBD[i] / M1);
                            g2 = (int)(ExBD[i] % M1);

                            d2f1 = (g1 >= EDGE && g1 < BoxShape[0] - EDGE) ? 
                                - KC1 * (F[(g1-2)*W1+g2] + F[(g1+2)*W1+g2])                                                                                       
                                + KC2 * (F[(g1-1)*W1+g2] + F[(g1+1)*W1+g2])
                                - KC3 * F[g1*W1+g2] : 
                                + F[(g1-1)*W1+g2] + F[(g1+1)*W1+g2] - 2.0 * F[g1*W1+g2];

                            d2f2 = (g2 >= EDGE && g2 < BoxShape[1] - EDGE) ? 
                                - KC1 * (F[g1*W1+(g2-2)] + F[g1*W1+(g2+2)])
                                + KC2 * (F[g1*W1+(g2-1)] + F[g1*W1+(g2+1)]) 
                                - KC3 * F[g1*W1+g2] : 
                                + F[g1*W1+(g2-1)] + F[g1*W1+(g2+1)] - 2.0 * F[g1*W1+g2];

                            KK1[g1*W1+g2] = I * (kh2m1hisq0 * d2f1 + kh2m2hisq1 * d2f2)  \
                                            - I * k2hb * POT[g1*W1+g2] * F[g1*W1+g2];

                            FF[g1*W1+g2] = F[g1*W1+g2] + KK1[g1*W1+g2] / 6.0;
                        }
                        #pragma omp single nowait
                        {
                            t_1_end = omp_get_wtime();
                            t_1_elapsed = t_1_end - t_1_begin;
                            t_overhead += t_1_elapsed;
                            if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-11: CASE 1 KK1) = %.4e sec\n", t_1_elapsed);
                            t_1_begin = omp_get_wtime();
                        }
                        // RK4-2
                        #pragma omp for private(g1, g2, d2f1, d2f2) schedule(static,4)
                        for (int i = 0; i < ExBD.size(); i++)  {

                            g1 = (int)(ExBD[i] / M1);
                            g2 = (int)(ExBD[i] % M1);

                            d2f1 = (g1 >= EDGE && g1 < BoxShape[0] - EDGE) ?
                                - KC1 * (F[(g1-2)*W1+g2] + F[(g1+2)*W1+g2] + 0.5 * (KK1[(g1-2)*W1+g2] + KK1[(g1+2)*W1+g2]))
                                + KC2 * (F[(g1-1)*W1+g2] + F[(g1+1)*W1+g2] + 0.5 * (KK1[(g1-1)*W1+g2] + KK1[(g1+1)*W1+g2]))
                                - KC3 * (F[g1*W1+g2] + 0.5 * KK1[g1*W1+g2]) :
                                + (F[(g1-1)*W1+g2] + F[(g1+1)*W1+g2] + 0.5 * (KK1[(g1-1)*W1+g2] + KK1[(g1+1)*W1+g2])) 
                                - 2.0 * F[g1*W1+g2] - KK1[g1*W1+g2];
                                
                            d2f2 = (g2 >= EDGE && g2 < BoxShape[1] - EDGE) ?
                                - KC1 * (F[g1*W1+(g2-2)] + F[g1*W1+(g2+2)] + 0.5 * (KK1[g1*W1+(g2-2)] + KK1[g1*W1+(g2+2)]))
                                + KC2 * (F[g1*W1+(g2-1)] + F[g1*W1+(g2+1)] + 0.5 * (KK1[g1*W1+(g2-1)] + KK1[g1*W1+(g2+1)]))
                                - KC3 * (F[g1*W1+g2] + 0.5 * KK1[g1*W1+g2]) :
                                + (F[g1*W1+(g2-1)] + F[g1*W1+(g2+1)] + 0.5 * (KK1[g1*W1+(g2-1)] + KK1[g1*W1+(g2+1)]))
                                - 2.0 * F[g1*W1+g2] - KK1[g1*W1+g2];

                            KK2[g1*W1+g2] = I * (kh2m1hisq0 * d2f1 + kh2m2hisq1 * d2f2)  \
                                            - I * k2hb * POT[g1*W1+g2] * (F[g1*W1+g2] + 0.5 * KK1[g1*W1+g2]);

                            FF[g1*W1+g2] += KK2[g1*W1+g2] / 3.0;
                        }
                        #pragma omp single nowait
                        {
                            t_1_end = omp_get_wtime();
                            t_1_elapsed = t_1_end - t_1_begin;
                            t_overhead += t_1_elapsed;
                            if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-12: CASE 1 KK2) = %.4e sec\n", t_1_elapsed);
                            t_1_begin = omp_get_wtime();
                        }
                        // RK4-3
                        #pragma omp for private(g1, g2, d2f1, d2f2) schedule(static,4)
                        for (int i = 0; i < ExBD.size(); i++)  {

                            g1 = (int)(ExBD[i] / M1);
                            g2 = (int)(ExBD[i] % M1);

                            d2f1 = (g1 >= EDGE && g1 < BoxShape[0] - EDGE) ?
                                - KC1 * (F[(g1-2)*W1+g2] + F[(g1+2)*W1+g2] + 0.5 * (KK2[(g1-2)*W1+g2] + KK2[(g1+2)*W1+g2]))
                                + KC2 * (F[(g1-1)*W1+g2] + F[(g1+1)*W1+g2] + 0.5 * (KK2[(g1-1)*W1+g2] + KK2[(g1+1)*W1+g2]))
                                - KC3 * (F[g1*W1+g2] + 0.5 * KK2[g1*W1+g2]) :
                                + (F[(g1-1)*W1+g2] + F[(g1+1)*W1+g2] + 0.5 * (KK2[(g1-1)*W1+g2] + KK2[(g1+1)*W1+g2]))
                                - 2.0 * F[g1*W1+g2] - KK2[g1*W1+g2];

                            d2f2 = (g2 >= EDGE && g2 < BoxShape[1] - EDGE) ?
                                - KC1 * (F[g1*W1+(g2-2)] + F[g1*W1+(g2+2)] + 0.5 * (KK2[g1*W1+(g2-2)] + KK2[g1*W1+(g2+2)]))
                                + KC2 * (F[g1*W1+(g2-1)] + F[g1*W1+(g2+1)] + 0.5 * (KK2[g1*W1+(g2-1)] + KK2[g1*W1+(g2+1)]))
                                - KC3 * (F[g1*W1+g2] + 0.5 * KK2[g1*W1+g2]) :
                                + (F[g1*W1+(g2-1)] + F[g1*W1+(g2+1)] + 0.5 * (KK2[g1*W1+(g2-1)] + KK2[g1*W1+(g2+1)]))
                                - 2.0 * F[g1*W1+g2] - KK2[g1*W1+g2];

                            KK3[g1*W1+g2] = I * (kh2m1hisq0 * d2f1 + kh2m2hisq1 * d2f2) \
                                            - I * k2hb * POT[g1*W1+g2] * (F[g1*W1+g2] + 0.5 * KK2[g1*W1+g2]);

                            FF[g1*W1+g2] += KK3[g1*W1+g2] / 3.0;
                        }
                        #pragma omp single nowait
                        {
                            t_1_end = omp_get_wtime();
                            t_1_elapsed = t_1_end - t_1_begin;
                            t_overhead += t_1_elapsed;
                            if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-13: CASE 1 KK3) = %.4e sec\n", t_1_elapsed);
                            t_1_begin = omp_get_wtime();
                        }
                        // RK4-4
                        #pragma omp for private(g1, g2, d2f1, d2f2) schedule(static,4)
                        for (int i = 0; i < ExBD.size(); i++)  {

                            g1 = (int)(ExBD[i] / M1);
                            g2 = (int)(ExBD[i] % M1);

                            d2f1 = (g1 >= EDGE && g1 < BoxShape[0] - EDGE) ?
                                - KC1 * (F[(g1-2)*W1+g2] + F[(g1+2)*W1+g2] + KK3[(g1-2)*W1+g2] + KK3[(g1+2)*W1+g2])
                                + KC2 * (F[(g1-1)*W1+g2] + F[(g1+1)*W1+g2] + KK3[(g1-1)*W1+g2] + KK3[(g1+1)*W1+g2])
                                - KC3 * (F[g1*W1+g2] + KK3[g1*W1+g2]) :
                                + (F[(g1-1)*W1+g2] + F[(g1+1)*W1+g2] + KK3[(g1-1)*W1+g2] + KK3[(g1+1)*W1+g2])
                                - 2.0 * (F[g1*W1+g2] + KK3[g1*W1+g2]);

                            d2f2 = (g2 >= EDGE && g2 < BoxShape[1] - EDGE) ?
                                - KC1 * (F[g1*W1+(g2-2)] + F[g1*W1+(g2+2)] + KK3[g1*W1+(g2-2)] + KK3[g1*W1+(g2+2)])
                                + KC2 * (F[g1*W1+(g2-1)] + F[g1*W1+(g2+1)] + KK3[g1*W1+(g2-1)] + KK3[g1*W1+(g2+1)])
                                - KC3 * (F[g1*W1+g2] + KK3[g1*W1+g2]) :
                                + (F[g1*W1+(g2-1)] + F[g1*W1+(g2+1)] + KK3[g1*W1+(g2-1)] + KK3[g1*W1+(g2+1)])
                                - 2.0 * (F[g1*W1+g2] + KK3[g1*W1+g2]);

                            KK4[g1*W1+g2] = I * (kh2m1hisq0 * d2f1 + kh2m2hisq1 * d2f2) \
                                            - I * k2hb * POT[g1*W1+g2] * (F[g1*W1+g2] + KK3[g1*W1+g2]);

                            FF[g1*W1+g2] += KK4[g1*W1+g2] / 6.0;
                        }
                        #pragma omp single nowait
                        {
                            t_1_end = omp_get_wtime();
                            t_1_elapsed = t_1_end - t_1_begin;
                            t_overhead += t_1_elapsed;
                            if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-14: CASE 1 KK4) = %.4e sec\n", t_1_elapsed);
                        }
                    }
                }  // if ( ExFF.size() > 0 )
            }  // if not ( isFirstExtrp ) and ( ExFF.size() == 1 )

            // Check Multiple Expanding 
            // TBL = index of FF that FF(TBL) is higher than TolL

            if ( ExFF.size() > 0 )  {

                t_1_begin = omp_get_wtime();

                TBL.clear();
                tmpVec.clear();

                #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2,f0,b1,b2,b3,b4,b5) schedule(static,4)
                for (int i = 0; i < ExFF.size(); i++)
                {
                    if (Check[i] == 1)  {
                        g1 = (int)(ExFF[i] / M1);
                        g2 = (int)(ExFF[i] % M1);

                        f0 = FF[g1*W1+g2];
                        b1 = std::norm(f0) >= TolH;
                        b2 = std::abs(FF[(g1+1)*W1+g2] - FF[(g1-1)*W1+g2]) >= TolHd2H0;
                        b3 = std::abs(FF[g1*W1+(g2+1)] - FF[g1*W1+(g2-1)]) >= TolHd2H1;
                        b4 = g1 > EDGE && g2 > EDGE;
                        b5 = (g1 < BoxShape[0]-EDGE-1) && (g2 < BoxShape[1]-EDGE-1);

                        if ((b1||b2||b3) && b4 && b5)  {
                            tmpVec.push_back(ExFF[i]);
                        }
                    }
                }       
                tmpVec.swap(TBL);
                tmpVec.clear(); 

                // TBL & TBL_P set difference
                tmpVec.resize(TBL_P.size() + TBL.size());
                QTR_SORT(TBL.begin(), TBL.end());
                QTR_SORT(TBL_P.begin(), TBL_P.end());
                it=std::set_difference( TBL.begin(), TBL.end(), TBL_P.begin(), TBL_P.end(), tmpVec.begin() );
                tmpVec.resize(it - tmpVec.begin()); 
                tmpVec.swap(TBL);
                tmpVec.clear();

                // Combine TBL and TBL_P
                TBL_P.reserve(TBL_P.size() + TBL.size());
                TBL_P.insert(TBL_P.end(), TBL.begin(), TBL.end());

                // Find unique elements
                QTR_SORT(TBL_P.begin(),TBL_P.end());
                it = std::unique (TBL_P.begin(), TBL_P.end()); 
                TBL_P.resize(std::distance(TBL_P.begin(),it));

                // Update isFirstExtrp
                Excount += 1;
                if (Excount == ExLimit) TBL.clear();

                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_overhead += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-c-3 CASE 1 TBL) = %.4e sec\n", t_1_elapsed); 
            }

        }  // while ( !isFullGrid && TBL.size() != 0  && Excount < ExLimit )
        
        // .........................................................................................

        // CASE 2: Truncating without extrapolation

        if ( !isExtrapolate && !isFullGrid )
        {
            #pragma omp parallel
            {
                #pragma omp single nowait
                {
                    t_1_begin = omp_get_wtime();
                }

                // RK4-1
                #pragma omp for private(d2f1, d2f2) schedule(static,4)
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        if (TAMask[i1*W1+i2])  {
                            d2f1 = (i1 >= EDGE && i1 < BoxShape[0] - EDGE) ? 
                            - KC1 * (F[(i1-2)*W1+i2] + F[(i1+2)*W1+i2])
                            + KC2 * (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2])
                            - KC3 * F[i1*W1+i2] : 
                            + F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] - 2.0 * F[i1*W1+i2];

                            d2f2 = (i2 >= EDGE && i2 < BoxShape[1] - EDGE) ? 
                            - KC1 * (F[i1*W1+(i2-2)] + F[i1*W1+(i2+2)])
                            + KC2 * (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)]) 
                            - KC3 * F[i1*W1+i2] : 
                            + F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] - 2.0 * F[i1*W1+i2];

                            KK1[i1*W1+i2] = I * (kh2m1hisq0 * d2f1 + kh2m2hisq1 * d2f2)  \
                                            - I * k2hb * POT[i1*W1+i2] * F[i1*W1+i2];

                            FF[i1*W1+i2] = F[i1*W1+i2] + KK1[i1*W1+i2] / 6.0;
                        }
                    }
                }

                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-21: CASE 2 KK1) = %.4e sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }

                // RK4-2
                #pragma omp for private(d2f1, d2f2) schedule(static,4)
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        if (TAMask[i1*W1+i2])  {
                            d2f1 = (i1 >= EDGE && i1 < BoxShape[0] - EDGE) ?
                                - KC1 * (F[(i1-2)*W1+i2] + F[(i1+2)*W1+i2] + 0.5 * (KK1[(i1-2)*W1+i2] + KK1[(i1+2)*W1+i2]))
                                + KC2 * (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] + 0.5 * (KK1[(i1-1)*W1+i2] + KK1[(i1+1)*W1+i2]))
                                - KC3 * (F[i1*W1+i2] + 0.5 * KK1[i1*W1+i2]) :
                                + (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] + 0.5 * (KK1[(i1-1)*W1+i2] + KK1[(i1+1)*W1+i2])) 
                                - 2.0 * F[i1*W1+i2] - KK1[i1*W1+i2];
                            
                            d2f2 = (i2 >= EDGE && i2 < BoxShape[1] - EDGE) ?
                                - KC1 * (F[i1*W1+(i2-2)] + F[i1*W1+(i2+2)] + 0.5 * (KK1[i1*W1+(i2-2)] + KK1[i1*W1+(i2+2)]))
                                + KC2 * (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] + 0.5 * (KK1[i1*W1+(i2-1)] + KK1[i1*W1+(i2+1)]))
                                - KC3 * (F[i1*W1+i2] + 0.5 * KK1[i1*W1+i2]) :
                                + (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] + 0.5 * (KK1[i1*W1+(i2-1)] + KK1[i1*W1+(i2+1)]))
                                - 2.0 * F[i1*W1+i2] - KK1[i1*W1+i2];

                            KK2[i1*W1+i2] = I * (kh2m1hisq0 * d2f1 + kh2m2hisq1 * d2f2) \
                                            - I * k2hb * POT[i1*W1+i2] * (F[i1*W1+i2] + 0.5 * KK1[i1*W1+i2]);

                            FF[i1*W1+i2] += KK2[i1*W1+i2] / 3.0;
                        }
                    }
                }

                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-22: CASE 2 KK2) = %.4e sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }

                // RK4-3
                #pragma omp for private(d2f1, d2f2) schedule(static,4)
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        if (TAMask[i1*W1+i2])  {
                            d2f1 = (i1 >= EDGE && i1 < BoxShape[0] - EDGE) ?
                                - KC1 * (F[(i1-2)*W1+i2] + F[(i1+2)*W1+i2] + 0.5 * (KK2[(i1-2)*W1+i2] + KK2[(i1+2)*W1+i2]))
                                + KC2 * (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] + 0.5 * (KK2[(i1-1)*W1+i2] + KK2[(i1+1)*W1+i2]))
                                - KC3 * (F[i1*W1+i2] + 0.5 * KK2[i1*W1+i2]) :
                                + (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] + 0.5 * (KK2[(i1-1)*W1+i2] + KK2[(i1+1)*W1+i2]))
                                - 2.0 * F[i1*W1+i2] - KK2[i1*W1+i2];

                            d2f2 = (i2 >= EDGE && i2 < BoxShape[1] - EDGE) ?
                                - KC1 * (F[i1*W1+(i2-2)] + F[i1*W1+(i2+2)] + 0.5 * (KK2[i1*W1+(i2-2)] + KK2[i1*W1+(i2+2)]))
                                + KC2 * (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] + 0.5 * (KK2[i1*W1+(i2-1)] + KK2[i1*W1+(i2+1)]))
                                - KC3 * (F[i1*W1+i2] + 0.5 * KK2[i1*W1+i2]) :
                                + (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] + 0.5 * (KK2[i1*W1+(i2-1)] + KK2[i1*W1+(i2+1)]))
                                - 2.0 * F[i1*W1+i2] - KK2[i1*W1+i2];

                            KK3[i1*W1+i2] = I * (kh2m1hisq0 * d2f1 + kh2m2hisq1 * d2f2) \
                                            - I * k2hb * POT[i1*W1+i2] * (F[i1*W1+i2] + 0.5 * KK2[i1*W1+i2]);

                            FF[i1*W1+i2] += KK3[i1*W1+i2] / 3.0;
                        }
                    }
                }

                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-23: CASE 2 KK3) = %.4e sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }

                // RK4-4
                #pragma omp for private(d2f1, d2f2) schedule(static,4)
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        if (TAMask[i1*W1+i2])  {
                            d2f1 = (i1 >= EDGE && i1 < BoxShape[0] - EDGE) ?
                                - KC1 * (F[(i1-2)*W1+i2] + F[(i1+2)*W1+i2] + KK3[(i1-2)*W1+i2] + KK3[(i1+2)*W1+i2])
                                + KC2 * (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] + KK3[(i1-1)*W1+i2] + KK3[(i1+1)*W1+i2])
                                - KC3 * (F[i1*W1+i2] + KK3[i1*W1+i2]) :
                                + (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] + KK3[(i1-1)*W1+i2] + KK3[(i1+1)*W1+i2])
                                - 2.0 * (F[i1*W1+i2] + KK3[i1*W1+i2]);

                            d2f2 = (i2 >= EDGE && i2 < BoxShape[1] - EDGE) ?
                                - KC1 * (F[i1*W1+(i2-2)] + F[i1*W1+(i2+2)] + KK3[i1*W1+(i2-2)] + KK3[i1*W1+(i2+2)])
                                + KC2 * (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] + KK3[i1*W1+(i2-1)] + KK3[i1*W1+(i2+1)])
                                - KC3 * (F[i1*W1+i2] + KK3[i1*W1+i2]) :
                                + (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] + KK3[i1*W1+(i2-1)] + KK3[i1*W1+(i2+1)])
                                - 2.0 * (F[i1*W1+i2] + KK3[i1*W1+i2]);

                            KK4[i1*W1+i2] = I * (kh2m1hisq0 * d2f1 + kh2m2hisq1 * d2f2) \
                                            - I * k2hb * POT[i1*W1+i2] * (F[i1*W1+i2] + KK3[i1*W1+i2]);

                            FF[i1*W1+i2] += KK4[i1*W1+i2] / 6.0;
                        }
                    }
                }

                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-24: CASE 2 KK4) = %.4e sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }
            } // OMP PARALLEL
        } 
        else if ( !isExtrapolate && isFullGrid )
        {
            // .........................................................................................

            // CASE 3: Full grid

            #pragma omp parallel
            {
                #pragma omp single nowait
                {
                    t_1_begin = omp_get_wtime();
                }
                // RK4-1
                #pragma omp for private(d2f1, d2f2) schedule(static,4)
                for (int i1 = EDGE; i1 < BoxShape[0] - EDGE; i1 ++)  {
                    for (int i2 = EDGE; i2 < BoxShape[1] - EDGE; i2 ++)  {
                        d2f1 = (i1 >= EDGE && i1 < BoxShape[0] - EDGE) ? 
                            - KC1 * (F[(i1-2)*W1+i2] + F[(i1+2)*W1+i2])
                            + KC2 * (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2])
                            - KC3 * F[i1*W1+i2] : 
                            + F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] - 2.0 * F[i1*W1+i2];

                        d2f2 = (i2 >= EDGE && i2 < BoxShape[1] - EDGE) ? 
                            - KC1 * (F[i1*W1+(i2-2)] + F[i1*W1+(i2+2)])
                            + KC2 * (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)]) 
                            - KC3 * F[i1*W1+i2] : 
                            + F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] - 2.0 * F[i1*W1+i2];

                        KK1[i1*W1+i2] = I * (kh2m1hisq0 * d2f1 + kh2m2hisq1 * d2f2)  \
                                        - I * k2hb * POT[i1*W1+i2] * F[i1*W1+i2];

                        FF[i1*W1+i2] = F[i1*W1+i2] + KK1[i1*W1+i2] / 6.0;
                    }
                }

                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_full += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-31: CASE 3 KK1) = %.4e sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }

                // RK4-2
                #pragma omp for private(d2f1, d2f2) schedule(static,4)
                for (int i1 = EDGE; i1 < BoxShape[0] - EDGE; i1 ++)  {
                    for (int i2 = EDGE; i2 < BoxShape[1] - EDGE; i2 ++)  {
                        d2f1 = (i1 >= EDGE && i1 < BoxShape[0]- EDGE) ?
                            - KC1 * (F[(i1-2)*W1+i2] + F[(i1+2)*W1+i2] + 0.5 * (KK1[(i1-2)*W1+i2] + KK1[(i1+2)*W1+i2]))
                            + KC2 * (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] + 0.5 * (KK1[(i1-1)*W1+i2] + KK1[(i1+1)*W1+i2]))
                            - KC3 * (F[i1*W1+i2] + 0.5 * KK1[i1*W1+i2]) :
                            + (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] + 0.5 * (KK1[(i1-1)*W1+i2] + KK1[(i1+1)*W1+i2])) 
                            - 2.0 * F[i1*W1+i2] - KK1[i1*W1+i2];
                        
                        d2f2 = (i2 >= EDGE && i2 < BoxShape[1] - EDGE) ?
                            - KC1 * (F[i1*W1+(i2-2)] + F[i1*W1+(i2+2)] + 0.5 * (KK1[i1*W1+(i2-2)] + KK1[i1*W1+(i2+2)]))
                            + KC2 * (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] + 0.5 * (KK1[i1*W1+(i2-1)] + KK1[i1*W1+(i2+1)]))
                            - KC3 * (F[i1*W1+i2] + 0.5 * KK1[i1*W1+i2]) :
                            + (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] + 0.5 * (KK1[i1*W1+(i2-1)] + KK1[i1*W1+(i2+1)]))
                            - 2.0 * F[i1*W1+i2] - KK1[i1*W1+i2];

                        KK2[i1*W1+i2] = I * (kh2m1hisq0 * d2f1 + kh2m2hisq1 * d2f2) \
                                        - I * k2hb * POT[i1*W1+i2] * (F[i1*W1+i2] + 0.5 * KK1[i1*W1+i2]);

                        FF[i1*W1+i2] += KK2[i1*W1+i2] / 3.0;
                    }
                }

                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_full += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-32: CASE 3 KK2) = %.4e sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }

                // RK4-3
                #pragma omp for private(d2f1, d2f2) schedule(static,4)
                for (int i1 = EDGE; i1 < BoxShape[0] - EDGE; i1 ++)  {
                    for (int i2 = EDGE; i2 < BoxShape[1] - EDGE; i2 ++)  {
                        d2f1 = (i1 >= EDGE && i1 < BoxShape[0] - EDGE) ?
                            - KC1 * (F[(i1-2)*W1+i2] + F[(i1+2)*W1+i2] + 0.5 * (KK2[(i1-2)*W1+i2] + KK2[(i1+2)*W1+i2]))
                            + KC2 * (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] + 0.5 * (KK2[(i1-1)*W1+i2] + KK2[(i1+1)*W1+i2]))
                            - KC3 * (F[i1*W1+i2] + 0.5 * KK2[i1*W1+i2]) :
                            + (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] + 0.5 * (KK2[(i1-1)*W1+i2] + KK2[(i1+1)*W1+i2]))
                            - 2.0 * F[i1*W1+i2] - KK2[i1*W1+i2];

                        d2f2 = (i2 >= EDGE && i2 < BoxShape[1] - EDGE) ?
                            - KC1 * (F[i1*W1+(i2-2)] + F[i1*W1+(i2+2)] + 0.5 * (KK2[i1*W1+(i2-2)] + KK2[i1*W1+(i2+2)]))
                            + KC2 * (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] + 0.5 * (KK2[i1*W1+(i2-1)] + KK2[i1*W1+(i2+1)]))
                            - KC3 * (F[i1*W1+i2] + 0.5 * KK2[i1*W1+i2]) :
                            + (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] + 0.5 * (KK2[i1*W1+(i2-1)] + KK2[i1*W1+(i2+1)]))
                            - 2.0 * F[i1*W1+i2] - KK2[i1*W1+i2];

                        KK3[i1*W1+i2] = I * (kh2m1hisq0 * d2f1 + kh2m2hisq1 * d2f2) \
                                        - I * k2hb * POT[i1*W1+i2] * (F[i1*W1+i2] + 0.5 * KK2[i1*W1+i2]);

                        FF[i1*W1+i2] += KK3[i1*W1+i2] / 3.0;
                    }
                }

                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_full += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-33: CASE 3 KK3) = %.4e sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }

                // RK4-4
                #pragma omp for private(d2f1, d2f2) schedule(static,4)
                for (int i1 = EDGE; i1 < BoxShape[0] - EDGE; i1 ++)  {
                    for (int i2 = EDGE; i2 < BoxShape[1] - EDGE; i2 ++)  {
                        d2f1 = (i1 >= EDGE && i1 < BoxShape[0] - EDGE) ?
                            - KC1 * (F[(i1-2)*W1+i2] + F[(i1+2)*W1+i2] + KK3[(i1-2)*W1+i2] + KK3[(i1+2)*W1+i2])
                            + KC2 * (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] + KK3[(i1-1)*W1+i2] + KK3[(i1+1)*W1+i2])
                            - KC3 * (F[i1*W1+i2] + KK3[i1*W1+i2]) :
                            + (F[(i1-1)*W1+i2] + F[(i1+1)*W1+i2] + KK3[(i1-1)*W1+i2] + KK3[(i1+1)*W1+i2])
                            - 2.0 * (F[i1*W1+i2] + KK3[i1*W1+i2]);

                        d2f2 = (i2 >= EDGE && i2 < BoxShape[1] - EDGE) ?
                            - KC1 * (F[i1*W1+(i2-2)] + F[i1*W1+(i2+2)] + KK3[i1*W1+(i2-2)] + KK3[i1*W1+(i2+2)])
                            + KC2 * (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] + KK3[i1*W1+(i2-1)] + KK3[i1*W1+(i2+1)])
                            - KC3 * (F[i1*W1+i2] + KK3[i1*W1+i2]) :
                            + (F[i1*W1+(i2-1)] + F[i1*W1+(i2+1)] + KK3[i1*W1+(i2-1)] + KK3[i1*W1+(i2+1)])
                            - 2.0 * (F[i1*W1+i2] + KK3[i1*W1+i2]);

                        KK4[i1*W1+i2] = I * (kh2m1hisq0 * d2f1 + kh2m2hisq1 * d2f2) \
                                        - I * k2hb * POT[i1*W1+i2] * (F[i1*W1+i2] + KK3[i1*W1+i2]);

                        FF[i1*W1+i2] += KK4[i1*W1+i2] / 6.0;
                    }
                }

                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_full += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-34: CASE 3 KK4) = %.4e sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }
            }
        }
        
        // .........................................................................................

        //  Renormalization to Theoretical Remains

        if (isReNorm)  {

            t_1_begin = omp_get_wtime();

            if (!isFullGrid)
            {    
                x2_flux_max = std::min(idx_des_2_max, x2_max);
                x2_flux_min = std::max(idx_des_2_min, x2_min);

                x1_flux_max = std::min(idx_dis_1_max, x1_max);
                x1_flux_min = std::max(idx_dis_1_min, x1_min);
            }    
            else 
            {    
                x2_flux_max = idx_des_2_max;
                x2_flux_min = idx_des_2_min;

                x1_flux_max = idx_dis_1_max;
                x1_flux_min = idx_dis_1_min;
            }    

            flux_des = 0.0;
            flux_dis = 0.0;

            for (int i2 = x2_flux_min; i2 < x2_flux_max; i2 ++)  {
                d1f1 = + D1 * (FF[(idx_des_1-2)*W1+i2] - FF[(idx_des_1+2)*W1+i2])
                       - D2 * (FF[(idx_des_1-1)*W1+i2] - FF[(idx_des_1+1)*W1+i2]);
                flux_des += h2m1h0 * imag(conj(FF[(idx_des_1)*W1+i2]) * d1f1);
            }
                                                                  
            for (int i1 = x1_flux_min; i1 < x1_flux_max; i1 ++)  {
                d1f2 = + D1 * (FF[i1*W1+(idx_dis_2-2)] - FF[i1*W1+(idx_dis_2+2)]) 
                       - D2 * (FF[i1*W1+(idx_dis_2-1)] - FF[i1*W1+(idx_dis_2+1)]);
                flux_dis += h2m2h1 * imag(conj(FF[i1*W1+(idx_dis_2)]) * d1f2);
            }

            flux_des *= H[1];
            flux_dis *= H[0];

            prob_remain -= (flux_des + flux_dis) * kk;

            t_1_end = omp_get_wtime();

            t_1_elapsed = t_1_end - t_1_begin;

            t_full += t_1_elapsed;
            t_truncate += t_1_elapsed;

            if (!QUIET && TIMING)
            {
                log->log(
                    "Elapsed time (omp-e-1-0 ReNorm) = %.4e sec\n",
                    t_1_elapsed
                );
            }
        }

        // NORMALIZATION AND TRUNCATION

        t_1_begin = omp_get_wtime();

        if (DEBUG)  {
            log->log("[DEBUG-TA-3A] [%d, %d][%d, %d]\n", x1_min, x1_max, x2_min, x2_max);
        }

        // Normalization

        norm = 0.0;

        if (!isFullGrid)  {
            #pragma omp parallel for reduction (+:norm) schedule(static,4)
            for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                int base_id = i1 * W1;
                for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                    if (TAMask[base_id+i2]) {
                        const auto &z = FF[base_id+i2];
                        norm += z.real()*z.real() + z.imag()*z.imag();
                    }
                }
            }
        }  
        else  {
            #pragma omp parallel for reduction (+:norm) schedule(static,4)
            for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
                int base_id = i1*W1;
                for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                    const auto &z = FF[base_id+i2];
                    norm += z.real()*z.real() + z.imag()*z.imag();
                }
            }
        }
        norm *= H[0] * H[1];

        if ( (tt + 1) % PERIOD == 0 )
            log->log("[Scatter2d] Time = %lf fs, Norm = %.16e\n", ( tt + 1 ) * kk * 0.024188843265864, norm);

        if (isReNorm)  {
            norm = sqrt(prob_remain/norm);
        }
        else {
            norm = sqrt(1.0/norm);
        }

        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        t_full += t_1_elapsed;
        t_truncate += t_1_elapsed;
        if (!QUIET && TIMING) log->log("Elapsed time (omp-e-1-1 Norm) = %.4e sec\n", t_1_elapsed);
        t_1_begin = omp_get_wtime();

        if (!isFullGrid)  {
            #pragma omp parallel for private(val) schedule(static,4)
            for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                    if (TAMask[i1*W1+i2])  {
                        val = norm * FF[i1*W1+i2];
                        FF[i1*W1+i2] = val;
                        F[i1*W1+i2] = val;
                        PF[i1*W1+i2] = std::norm(F[i1*W1+i2]);
                    }
                }
            }
        }
        else  {
            #pragma omp parallel for private(val) schedule(static,4)
            for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
                for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                    val = norm * FF[i1*W1+i2];
                    FF[i1*W1+i2] = val;
                    F[i1*W1+i2] = val;
                    PF[i1*W1+i2] = std::norm(F[i1*W1+i2]);
                }
            }
        }

        /* -------------------------------------------- */

        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        t_full += t_1_elapsed;
        t_truncate += t_1_elapsed;
        if (!QUIET && TIMING) log->log("Elapsed time (omp-e-1-2 FF) = %.4e sec\n", t_1_elapsed); 

        // Truncation and TA

        if ( !isFullGrid )
        {
            t_1_begin = omp_get_wtime();

            //#pragma omp parallel for schedule(auto)
            #pragma omp parallel for schedule(static,4)
            for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                    if (TAMask[i1*W1+i2])  {
                        if (PF[i1*W1+i2] < TolH)  {
                            if (std::abs(FF[(i1+1)*W1+i2] - FF[(i1-1)*W1+i2]) < TolHd2H0) {
                                if (std::abs(FF[i1*W1+(i2+1)] - FF[i1*W1+(i2-1)]) < TolHd2H1) {
                                    F[i1*W1+i2] = xZERO;
                                }
                            }
                        }
                    }
                }
            } 
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-3-1 TA) = %.4e sec\n", t_1_elapsed);
            t_1_begin = omp_get_wtime();

            #pragma omp parallel for schedule(static,4)
            for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                    if (F[i1*W1+i2] == 0.0)  {
                        TAMask[i1*W1+i2] = 0;
                        PF[i1*W1+i2] = 0;
                    }  else  {
                        if (!TAMask[i1*W1+i2]) {
                            TAMask[i1*W1+i2] = 1;
                        }
                    }
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-3-2 TA) = %.4e sec\n", t_1_elapsed);

            // Rebuild TA box

            t_1_begin = omp_get_wtime();

            x1_min_tmp = (x1_min - ExLimit) < EDGE ? EDGE : (x1_min - ExLimit);
            x1_max_tmp = (x1_max + ExLimit) > BoxShape[0]-EDGE-1 ? BoxShape[0]-EDGE-1 : (x1_max + ExLimit);
            x2_min_tmp = (x2_min - ExLimit) < EDGE ? EDGE : (x2_min - ExLimit);
            x2_max_tmp = (x2_max + ExLimit) > BoxShape[1]-EDGE-1 ? BoxShape[1]-EDGE-1 : (x2_max + ExLimit);

            ta_size = 0;

            #pragma omp parallel for reduction(min: x1_min,x2_min) \
                                     reduction(max: x1_max,x2_max) \
                                     reduction(+: ta_size) schedule(static,4)
            for (int i1 = x1_min_tmp; i1 <= x1_max_tmp; i1 ++)  {
                for (int i2 = x2_min_tmp; i2 <= x2_max_tmp; i2 ++)  {
                    if (TAMask[i1*W1+i2])  {
                        if (i1 < x1_min)  x1_min = i1;
                        if (i1 > x1_max)  x1_max = i1;
                        if (i2 < x2_min)  x2_min = i2;
                        if (i2 > x2_max)  x2_max = i2;
                        ta_size += 1;
                    }
                }
            }

            if (DEBUG)
                log->log("[DEBUG-TA-%d] ta_size = %d\n",tt+1,ta_size);

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-4 TA rebuild) = %.4e sec\n", t_1_elapsed);

            // TB
            
            t_1_begin = omp_get_wtime();

            if (ta_size == 0)
                tb_size = 0;
            else  {
                #pragma omp parallel for reduction(merge: tmpVec) schedule(static,4)
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        if (TAMask[i1*W1+i2])  {
                            if (!TAMask[(i1-1)*W1+i2] || !TAMask[(i1+1)*W1+i2] || \
                                !TAMask[i1*W1+(i2-1)] || !TAMask[i1*W1+(i2+1)])
                                tmpVec.push_back(i1*W1+i2);
                        }
                    }
                }
                tmpVec.swap(TB);
                tmpVec.clear();
                tb_size = TB.size();
            }

            if (DEBUG)
                log->log("[DEBUG-TB-%d] tb_size = %d\n",tt+1,tb_size);

            /* -------------------------------------------- */

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-5 TB) = %.4e sec\n", t_1_elapsed);

            // TA expansion

            t_1_begin = omp_get_wtime();
            
            #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2) schedule(static,4)
            for (int i = 0; i < TB.size(); i++)
            {
                g1 = (int)(TB[i] / M1);
                g2 = (int)(TB[i] % M1);

                if ( g1+1 < BoxShape[0]-EDGE-1 && !TAMask[(g1+1)*W1+g2])
                    tmpVec.push_back(GridToIdx(g1+1,g2));
                if ( g1-1 > EDGE && !TAMask[(g1-1)*W1+g2])
                    tmpVec.push_back(GridToIdx(g1-1,g2));
                if ( g2+1 < BoxShape[1]-EDGE-1 && !TAMask[g1*W1+(g2+1)])
                    tmpVec.push_back(GridToIdx(g1,g2+1));
                if ( g2-1 > EDGE && !TAMask[g1*W1+(g2-1)])
                    tmpVec.push_back(GridToIdx(g1,g2-1));
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-6 TAEX-A) = %.4e sec\n", t_1_elapsed);
            t_1_begin = omp_get_wtime();

            #pragma omp parallel for reduction(min: x1_min,x2_min) \
                                     reduction(max: x1_max,x2_max) \
                                     private(g1,g2) schedule(static,4)
            for (int i = 0; i < tmpVec.size(); i ++)  {
                g1 = (int)(tmpVec[i] / M1);
                g2 = (int)(tmpVec[i] % M1);
                TAMask[tmpVec[i]] = 1;
                x1_min = (g1 < x1_min) ? g1 : x1_min;
                x1_max = (g1 > x1_max) ? g1 : x1_max;
                x2_min = (g2 < x2_min) ? g2 : x2_min;
                x2_max = (g2 > x2_max) ? g2 : x2_max;
            }
            tmpVec.clear();

            ta_size = 0;
            #pragma omp parallel for reduction(+: ta_size) schedule(static,4)
            for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                    if (TAMask[i1*W1+i2]) {
                        ta_size += 1; 
                    }
                }
            }

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-7 TARB-A) = %.4e sec\n", t_1_elapsed);

            if (DEBUG)  {
                log->log("[DEBUG-TA-4A] [%d, %d][%d, %d]\n", x1_min, x1_max, x2_min, x2_max);
            }

            if (isToggle && ((ta_size * 1.0 / GRIDS_TOT) > toggle_threshold ))  {
                isFullGrid = true;
                log->log("[Scatter2d] threshold = %lf, ta_A/ntotal = %lf\n", toggle_threshold, ta_size*1.0/GRIDS_TOT);
                log->log("[Scatter2d] Switching to full-grid mode at step: %d\n", tt + 1);
            }
        }

        if ( (tt + 1) % PERIOD == 0 )
        {   
            write_original_observables(tt + 1);
            t_0_end = omp_get_wtime();
            t_0_elapsed = t_0_end - t_0_begin;

            if ( !QUIET ) log->log("[Scatter2d] Step: %d, Elapsed time: %.4e sec\n", tt + 1, t_0_elapsed);

            if ( !isFullGrid && !QUIET )  {
                log->log("[Scatter2d] norm = %.2e\n", norm);
                log->log("[Scatter2d] TA size = %d, TB size = %d\n", ta_size, tb_size);
                log->log("[Scatter2d] TA Range [%d, %d][%d, %d]\n", x1_min, x1_max, x2_min, x2_max);
                log->log("[Scatter2d] TA / total grids = %lf\n", ( ta_size * 1.0 ) / GRIDS_TOT);
                log->log("[Scatter2d] ExCount = %d , ExLimit =  %d\n", Excount, ExLimit);
                log->log("[Scatter2d] Core computation time = %lf\n", t_truncate);
                log->log("[Scatter2d] Overhead time = %lf\n\n", t_overhead);
            }
            else if ( isFullGrid && !QUIET )  {

                log->log("[Scatter2d] Core computation time = %lf\n", t_full);
            }
            if ( !QUIET ) log->log("\n........................................................\n\n");
        }         
    } // Time iteration 

    if (flux_d_file) fclose(flux_d_file);
    if (flux_1_file) fclose(flux_1_file);
    if (flux_2_file) fclose(flux_2_file);
    if (flux_3_file) fclose(flux_3_file);
    if (prob_r_a_file) fclose(prob_r_a_file);
    if (prob_a_d_file) fclose(prob_a_d_file);
    if (prob_d_1_file) fclose(prob_d_1_file);
    if (prob_1_2_file) fclose(prob_1_2_file);
    if (prob_2_3_file) fclose(prob_2_3_file);
    if (prob_3_p_file) fclose(prob_3_p_file);

    delete[] F;
    delete[] FF;
    delete[] PF;
    delete[] POT;
    delete[] KK1;
    delete[] KK2;
    delete[] KK3;
    delete[] KK4;

    if (TAMask)
        delete[] TAMask;

    log->log("[Scatter2d] Evolve done.\n");
}
/* ------------------------------------------------------------------------------- */

inline int Scatter2d::GridToIdx(int x1, int x2)
{
    return (int)(x1 * W1 + x2);
}
/* ------------------------------------------------------------------------------- */

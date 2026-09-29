// ==============================================================================
//
//  Parameters.cpp
//  QTR
//
//  Created by Albert Lu on 8/4/18.
//  alu@tacc.utexas.edu
//
//  Last modified on 1/27/19
//
//  Note:
//
// ==============================================================================

#include <cstring>
#include <errno.h>
#include <float.h>
#include <time.h>

#include "Constants.h"
#include "InputFile.h"
#include "Log.h"
#include "Parameters.h"

using namespace QTR_NS;
using std::string;

/* ------------------------------------------------------------------------------- */

Parameters::Parameters(class QTR *q) : Pointers(q)
{
    log = qtr->log;
    
    // Set default values here
    
    // MAIN //
    job               = "point";
    inFilename        = "test.in";
    logFilename       = "test.log";
    quiet             = false;
    writeLog          = true;
    
    // MPI //
    me                = 0;
    i                 = 0;
    
    // RANDOM //
    rngSeed           = -1;
    rngType           = "mars";
}
/* ------------------------------------------------------------------------------- */

Parameters::~Parameters()
{
    return;
}
/* ------------------------------------------------------------------------------- */

string Parameters::toLowerCase(string s)
{
    for (string::size_type i = 0; i < s.length(); ++i)
    {
        s[i] = tolower(s[i]);
    }
    return s;
}
/* ------------------------------------------------------------------------------- */

int Parameters::load(string filename)
{
    FILE *fh;
    
    fh = fopen(filename.c_str(), "rb");
    
    if (fh == NULL) {
        fprintf(stderr, "error: %s\n", strerror(errno));
        return 1;
    }
    
    int error = load(fh);
    
    fclose(fh);
    
    return error;
}
/* ------------------------------------------------------------------------------- */

int Parameters::load(FILE *file){
        
    InputFile ini;
    ini.CaseInsensitive();
    
    int error = 0;
    
    if(ini.ReadFile(file))
    {
        // MAIN //
        job         = toLowerCase(ini.GetValue("MAIN", "job"));
        inFilename  = ini.GetValue ("MAIN", "in_filename" , inFilename);
        logFilename = ini.GetValue ("MAIN", "log_filename", logFilename);
        quiet       = ini.GetValueB("MAIN", "quiet", quiet);
        timing      = ini.GetValueB("MAIN", "timing", timing);
        writeLog    = ini.GetValueB("MAIN", "write_log", writeLog);
        // SCATTERXD //
        scxd_isDebug    = ini.GetValueB("SCATTERXD", "isDebug", 0);
        scxd_isFullGrid = ini.GetValueB("SCATTERXD", "isFullGrid", 1);  
        scxd_isTrans    = ini.GetValueB("SCATTERXD", "isTrans", 1);
        scxd_isAcf      = ini.GetValueB("SCATTERXD", "isAcf", 1);
        scxd_isFlux      = ini.GetValueB("SCATTERXD", "isFlux", 1);
        scxd_isCondProb      = ini.GetValueB("SCATTERXD", "isCondProb", 1);
        scxd_isLShape    = ini.GetValueB("SCATTERXD", "isLShape", 0);
        // Keep the original package's norm-preserving behavior by default.
        // The theoretical-remain re-normalization is an opt-in target feature.
        scxd_isReNorm      = ini.GetValueB("SCATTERXD", "isReNorm", 0);
        scxd_isFluxNProb      = ini.GetValueB("SCATTERXD", "isFluxNProb", 1);
        scxd_isPrintEdge = ini.GetValueB("SCATTERXD", "isPrintEdge", 0);
        scxd_isPrintDensity = ini.GetValueB("SCATTERXD", "isPrintDensity", 0);
        scxd_isToggle        = ini.GetValueB("SCATTERXD", "isToggle", 0);
        scxd_isPrintWavefunc = ini.GetValueB("SCATTERXD", "isPrintWavefunc", 0);
        scxd_dimensions = ini.GetValueI("SCATTERXD", "dimensions", 2);  
        scxd_period = ini.GetValueI("SCATTERXD", "period", 100);
        scxd_sortperiod = ini.GetValueI("SCATTERXD", "sortperiod", 100);
        scxd_printperiod = ini.GetValueI("SCATTERXD", "printperiod", 100);
        scxd_k      = ini.GetValueF("SCATTERXD", "k", 0.001);
        scxd_h1     = ini.GetValueF("SCATTERXD", "h1", 0.1);
        scxd_h2     = ini.GetValueF("SCATTERXD", "h2", 0.1);
        scxd_h3     = ini.GetValueF("SCATTERXD", "h3", 0.1);
        scxd_h4     = ini.GetValueF("SCATTERXD", "h4", 0.1);
        scxd_xi1    = ini.GetValueF("SCATTERXD", "xi1", 0.1);
        scxd_xf1    = ini.GetValueF("SCATTERXD", "xf1", 0.1);
        scxd_xi2    = ini.GetValueF("SCATTERXD", "xi2", 0.1);
        scxd_xf2    = ini.GetValueF("SCATTERXD", "xf2", 0.1);
        scxd_xi3    = ini.GetValueF("SCATTERXD", "xi3", 0.1);
        scxd_xf3    = ini.GetValueF("SCATTERXD", "xf3", 0.1);
        scxd_xi4    = ini.GetValueF("SCATTERXD", "xi4", 0.1);
        scxd_xf4    = ini.GetValueF("SCATTERXD", "xf4", 0.1);
        scxd_l_shape_x1_cut = ini.GetValueF("SCATTERXD", "l_shape_x1_cut", scxd_xf1);
        scxd_l_shape_x2_cut = ini.GetValueF("SCATTERXD", "l_shape_x2_cut", scxd_xf2);
        scxd_Tf     = ini.GetValueF("SCATTERXD", "Tf", 0.1);
        scxd_x01    = ini.GetValueF("SCATTERXD", "x01", 0);
        scxd_x02    = ini.GetValueF("SCATTERXD", "x02", 0);
        scxd_x03    = ini.GetValueF("SCATTERXD", "x03", 0);
        scxd_x04    = ini.GetValueF("SCATTERXD", "x04", 0);
        scxd_a1     = ini.GetValueF("SCATTERXD", "a1", 0);
        scxd_a2     = ini.GetValueF("SCATTERXD", "a2", 0);
        scxd_a3     = ini.GetValueF("SCATTERXD", "a3", 0);
        scxd_a4     = ini.GetValueF("SCATTERXD", "a4", 0);
        scxd_p1     = ini.GetValueF("SCATTERXD", "p1", 0);
        scxd_p2     = ini.GetValueF("SCATTERXD", "p2", 0);
        scxd_p3     = ini.GetValueF("SCATTERXD", "p3", 0);  
        scxd_p4     = ini.GetValueF("SCATTERXD", "p4", 0);                
        scxd_hb     = ini.GetValueF("SCATTERXD", "hd", 1);
        scxd_m      = ini.GetValueF("SCATTERXD", "m",  1);
        scxd_m1      = ini.GetValueF("SCATTERXD", "m1",  1.0);
        scxd_m2      = ini.GetValueF("SCATTERXD", "m2",  1.0);
        scxd_TolH     = ini.GetValueF("SCATTERXD", "TolH", 0);
        scxd_TolL     = ini.GetValueF("SCATTERXD", "TolL", 0);
        scxd_TolHd    = ini.GetValueF("SCATTERXD", "TolHd", 0);
        scxd_TolLd    = ini.GetValueF("SCATTERXD", "TolLd", 0);
        scxd_ExReduce = ini.GetValueF("SCATTERXD", "ExReduce", 0);
	    scxd_ExLimit = ini.GetValueI("SCATTERXD", "ExLimit", 3);
        scxd_Vmode_1  = ini.GetValueI("SCATTERXD", "Vmode_1", 0);
        scxd_Vmode_2  = ini.GetValueI("SCATTERXD", "Vmode_2", 0);
        scxd_Vmode_3  = ini.GetValueI("SCATTERXD", "Vmode_3", 0);
        scxd_Vmode_4  = ini.GetValueI("SCATTERXD", "Vmode_4", 0);
        scxd_w        = ini.GetValueF("SCATTERXD", "w", 1.0);    // HO specific
        scxd_V0       = ini.GetValueF("SCATTERXD", "V0",   1.0); // Eckart potential 
        scxd_ek2v     = ini.GetValueF("SCATTERXD", "ek2v", 1.0); 
        scxd_alpha    = ini.GetValueF("SCATTERXD", "Alpha", 0.4);           
        scxd_k0     = ini.GetValueF("SCATTERXD", "k0",  0.09);   // Related HO
        scxd_sig    = ini.GetValueF("SCATTERXD", "sig", 0.1);        
        scxd_lan    = ini.GetValueF("SCATTERXD", "lan", 1.0);          
        scxd_De     = ini.GetValueF("SCATTERXD", "De", 40.0);    // Morse
        scxd_Da     = ini.GetValueF("SCATTERXD", "Da", 0.25);        
        scxd_r0     = ini.GetValueF("SCATTERXD", "r0", 0.0); 
        scxd_lambda = ini.GetValueF("SCATTERXD", "lambda", 0.11803);  // Henon-Heiles
        scxd_sigma  = ini.GetValueF("SCATTERXD", "sigma", 1.0);       // Seck
        scxd_beta   = ini.GetValueF("SCATTERXD", "beta", 1.0);        // Gauk   
        scxd_dk     = ini.GetValueF("SCATTERXD", "dk", 0.01);         // Spectrum dk
        scxd_kmax   = ini.GetValueF("SCATTERXD", "kmax", 10.0);       // Spectrum kmax
        scxd_idx_abs_1 = ini.GetValueI("SCATTERXD", "idx_abs_1", 100);
        scxd_idx_dsc_1 = ini.GetValueI("SCATTERXD", "idx_dsc_1", 54);
        scxd_idx_dsc_2 = ini.GetValueI("SCATTERXD", "idx_dsc_2", 52);
        scxd_idx_df1_1 = ini.GetValueI("SCATTERXD", "idx_df1_1", 60);
        scxd_idx_df1_2 = ini.GetValueI("SCATTERXD", "idx_df1_2", 205);
        scxd_idx_df2_1 = ini.GetValueI("SCATTERXD", "idx_df2_1", 60);
        scxd_idx_df2_2 = ini.GetValueI("SCATTERXD", "idx_df2_2", 409);
        scxd_idx_df3_1 = ini.GetValueI("SCATTERXD", "idx_df3_1", 60);
        scxd_idx_df3_2 = ini.GetValueI("SCATTERXD", "idx_df3_2", 609);
        scxd_idx_des_2_min = ini.GetValueI("SCATTERXD", "idx_des_2_min", 0); 
        scxd_idx_des_2_max = ini.GetValueI("SCATTERXD", "idx_des_2_max", 330); 
        scxd_idx_des_1 = ini.GetValueI("SCATTERXD", "idx_des_1", 451); 
        scxd_idx_dis_1_min = ini.GetValueI("SCATTERXD", "idx_dis_1_min", 0); 
        scxd_idx_dis_1_max = ini.GetValueI("SCATTERXD", "idx_dis_1_max", 550); 
        scxd_idx_dis_2 = ini.GetValueI("SCATTERXD", "idx_dis_2", 285); 
        scxd_edge   = ini.GetValueI("SCATTERXD", "edge", 2);          // Edge size
        scxd_toggle_threshold = ini.GetValueF("SCATTERXD", "toggle_threshold", 1.0); 
       
        // RANDOM //
        rngSeed     = ini.GetValueL("RANDOM", "random_seed" , rngSeed);
        rngType     = ini.GetValue ("RANDOM", "random_type" , rngType);
    }
    else
    {
        fprintf(stderr, "Couldn't parse the ini file.\n");
        error = 1;
    }
    
    return error;
}
/* ------------------------------------------------------------------------------- */

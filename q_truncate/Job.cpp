// ==============================================================================
//
//  Job.cpp
//  QTR
//
//  Created by Albert Lu on 11/2/20.
//  alu@tacc.utexas.edu
//
//  Last modified on 4/30/26
//
//  Note:
//
// ==============================================================================

# include "Job.h"
# include "Parameters.h"

# include "Job_Scatter2d.h"

using namespace QTR_NS;

const char Job::SCATTER2D[] = "scatter2d";

/* ------------------------------------------------------------------------------- */

Job *Job::getJob(class QTR *qtr) {
    
    Job *job = NULL;
    
    if (qtr->parameters->job == SCATTER2D)
    {
        job = new JobScatter2d(qtr);
    }
    return job;
}
/* ----------------------------------------------------------------- */

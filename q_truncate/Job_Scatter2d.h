// ==============================================================================
//
//  Job_Scatter2d.h
//  QTR
//
//  Created by Albert Lu on 1/14/21.
//  alu@tacc.utexas.edu
//
//  Last modified on 4/30/26
//
//  Note:
//
// ==============================================================================

#ifndef QTR_JOB_SCATTER2D_H
#define QTR_JOB_SCATTER2D_H

#include "Job.h"
#include "Scatter2d.h"

namespace QTR_NS  {

    class JobScatter2d: public Job    {
        
    public:
        JobScatter2d(class QTR *);
        ~JobScatter2d(void);
        
        void run(class QTR *);
        
    private:
        
        Log *log;
        Parameters *parameters;
        QTR *qtr;
        Scatter2d  *scat;
    };
}
#endif /* QTR_JOB_SCATTER2D_H */

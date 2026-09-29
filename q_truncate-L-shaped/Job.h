// ==============================================================================
//
//  Job.h
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

#ifndef QTR_JOB_H
#define QTR_JOB_H

#include "Qtr.h"

namespace QTR_NS {
    
    class Job {
        
    public:
        
        virtual ~Job() {};
        
        virtual void run(class QTR *) = 0;
        
        static Job *getJob(class QTR *);

        static const char SCATTER2D[];
    };
}
#endif /* QTR_JOB_H */

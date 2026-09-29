// ==============================================================================
//
//  Log.h
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

#ifndef QTR_LOG_H
#define QTR_LOG_H

#include <stdio.h>
#include <stdarg.h>

#include "Qtr.h"
#include "Pointers.h"

namespace QTR_NS {

    class Log : protected Pointers  {
        
    public:
        
        Log(class QTR *);
        ~Log();
        
        void log_init(class QTR *, char *filename);
        void log_close();
        void log(const char* format, ...);
        void log_file(const char* format, ...);
    };
    
}

#endif /* QTR_LOG_H */

// ==============================================================================
//
//  Error.h
//  QTR
//
//  Created by Albert Lu on 8/4/18.
//  alu@tacc.utexas.edu
//
//  Last modified on 10/31/18
//
//  Note:
//
// ==============================================================================

#ifndef QTR_ERROR_H
#define QTR_ERROR_H

#ifdef QTRMPI
#include <mpi.h>
#endif

#include <string>

#include "Pointers.h"

using std::string;

namespace QTR_NS   {
    
    class Error : protected Pointers   {
        
    public:  
        Error(class QTR *);
        ~Error();
        void abort_all();
        
# ifdef QTRMPI
        long     i;
        int      me;
        int      nprocs;
# endif
        
    private:  
        Log  *log;       
        char Error_1;
    };   
}

#endif /* QTR_ERROR_H */

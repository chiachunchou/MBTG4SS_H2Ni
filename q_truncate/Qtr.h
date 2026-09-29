// ==============================================================================
//
//  Qtr.h
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


#ifndef QTR_QTR_H
#define QTR_QTR_H

#ifdef QTRMPI
#include <mpi.h>
#endif

#include <stdio.h>

namespace QTR_NS  {

    class QTR {
        
    public:
        
# ifdef QTRMPI
        QTR(MPI_Comm communicator, char **arg);
        MPI_Comm world;
# else
        QTR(char **arg);
# endif
        
        ~QTR();
        
        void run();
    
        class Error      *error;
        class Job        *job;
        class Log        *log;
        class Parameters     *parameters;
        class RandNum      *randnum;
        
    private:
        char *inFilePtr;
    };
    
}

#endif /* QTR_QTR_H */

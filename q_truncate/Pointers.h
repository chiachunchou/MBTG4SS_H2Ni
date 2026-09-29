// ==============================================================================
//
//  Pointers.h
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

#ifndef QTR_POINTERS_H
#define QTR_POINTERS_H

#include "Qtr.h"

namespace QTR_NS {

    class Pointers {
        
    public:
        
        Pointers(QTR *ptr):
        qtr(ptr),
        error(ptr->error),
        log(ptr->log),
        parameters(ptr->parameters),
        randnum(ptr->randnum) {}
        
        virtual ~Pointers() {}
        
    protected:
        
        QTR    *qtr;
        Error    *&error;
        Log    *&log;
        Parameters *&parameters;
        RandNum  *&randnum;
    };

}

#endif /* QTR_POINTERS_H */

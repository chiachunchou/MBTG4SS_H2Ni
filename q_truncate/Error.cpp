// ==============================================================================
//
//  Error.cpp
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

#include "Error.h"
#include "Log.h"
#include "Parameters.h"

#include <iostream>

using std::cout;
using std::endl;

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

Error::Error(QTR *qtr) : Pointers(qtr)
{
    log = qtr->log;
    
    return;
}
/* ------------------------------------------------------------------------------- */

Error::~Error()  {
    
    return;
}
/* ------------------------------------------------------------------------------- */

void Error::abort_all()
{
    
#ifdef QTRMPI
    log->log("[Error] QTR aborting ...");
    MPI_Abort(qtr->parameters->world, 1);
#else
    log->log("[Error] QTR aborting ...");
    abort();
#endif
}
/* ------------------------------------------------------------------------------- */


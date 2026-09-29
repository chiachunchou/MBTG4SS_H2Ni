// ==============================================================================
//
//  Qtr.cpp
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

#ifdef QTRMPI
#include <mpi.h>
#endif

#include <string>

#include "Qtr.h"
#include "Error.h"
#include "Job.h"
#include "Log.h"
#include "Parameters.h"
#include "RandNum.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------- */

#ifdef QTRMPI

QTR::QTR(MPI_Comm communicator,char **arg)
    : error(nullptr), job(nullptr), log(nullptr), parameters(nullptr), randnum(nullptr), inFilePtr(arg[1])
{
    int     r,n;
    char    name[100];
    
    MPI_Group world_group, new_group;
    MPI_Comm  new_comm;
    
    MPI_Comm_rank(communicator,&r);
    MPI_Comm_size(communicator,&n);
        
    parameters = new Parameters(this);
    
    world = communicator;
    inFilePtr = arg[1];
    
    parameters->universe = communicator;
    parameters->me = r;
    parameters->nprocs = n;
    parameters->inFilename = std::string(arg[1]);
    
    sprintf(name,"log%d.txt",r);
    parameters->logFilename = std::string(name);
    
    log = new Log(this);
    error = new Error(this);
    
    MPI_Comm_group(communicator, &world_group);
    MPI_Group_incl(world_group, 1, &r, &new_group);
    MPI_Comm_create(communicator, new_group, &new_comm);
    
    if (new_comm != MPI_COMM_NULL)
    {
        parameters->world = new_comm;
    }
}
# else

/* ------------------------------------------------------------------------- */

QTR::QTR(char **arg)
    : error(nullptr), job(nullptr), log(nullptr), parameters(nullptr), randnum(nullptr), inFilePtr(arg[1])
{
    parameters = new Parameters(this);
    log = new Log(this);
    error = new Error(this);
}
# endif
/* ------------------------------------------------------------------------- */

QTR::~QTR()
{
    delete error;
    delete log;
    delete parameters;
}
/* ------------------------------------------------------------------------- */

void QTR::run()
{
    /* Load initial parameters */

    FILE *inFile;
    inFile = fopen(inFilePtr,"r");
    if (inFile == nullptr) {
        fprintf(stderr, "error: unable to open input file %s\n", inFilePtr);
        return;
    }
    if (this->parameters->load(inFile) != 0) {
        fclose(inFile);
        fprintf(stderr, "error: unable to parse input file %s\n", inFilePtr);
        return;
    }
    fclose(inFile);

    /* Activate Job */
    Job *active_job = Job::getJob(this);
    if (active_job == nullptr) {
        fprintf(stderr, "error: unsupported job '%s'\n", this->parameters->job.c_str());
        return;
    }
    active_job->run(this);
    delete active_job;
}
/* ------------------------------------------------------------------------- */

// ==============================================================================
//
//  Log.cpp
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

#include <stdio.h>
#include <stdarg.h>

#include "Log.h"
#include "Parameters.h"

using namespace QTR_NS;

static FILE *logfile = NULL;

/* ------------------------------------------------------------------------------- */

Log::Log(QTR *qtr) : Pointers(qtr)
{
    log_init(qtr, (char*)qtr->parameters->logFilename.c_str());
}
/* ------------------------------------------------------------------------------- */

Log::~Log()  {
    
    return;
}
/* ------------------------------------------------------------------------------- */

void Log::log_init(QTR *p, char *filename)
{
    qtr = p;
    
    if (logfile == NULL)
    {
        logfile = fopen(filename, "w");
    }
}
/* ------------------------------------------------------------------------------- */

void Log::log_close()
{
    if (logfile != NULL)
    {
        fclose(logfile);
        logfile = NULL;
    }
}
/* ------------------------------------------------------------------------------- */

void Log::log(const char* format, ...)
{
    if (logfile == NULL)  {
        fprintf(stderr, "error: log() called before log_init\n");
        return;
    }
    
    va_list args;
    
    if (qtr->parameters->quiet == false)   {
        va_start(args, format);
        vfprintf(stdout, format, args);
        va_end(args);
        //uncomment to ensure that info is saved to log immediately
        fflush(logfile);
    }
    
    if (qtr->parameters->writeLog == true)   {
        va_start(args, format);
        vfprintf(logfile, format, args);
        va_end(args);
        //uncomment to ensure that info is saved to log immediately
        fflush(logfile);
    }
}
/* ------------------------------------------------------------------------------- */

void Log::log_file(const char* format, ...)
{
    if (logfile == NULL)  {
        fprintf(stderr, "error: log() called before log_init\n");
        return;
    }
    
    if (qtr->parameters->writeLog == true)   {
        va_list args;
        va_start(args, format);
        vfprintf(logfile, format, args);
        va_end(args);
        //uncomment to ensure that info is saved to log immediately
        fflush(logfile);
    }
}
/* ------------------------------------------------------------------------------- */



// ==============================================================================
//
//  main.cpp
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

# include <cstdlib>
# include <iostream>

# ifdef QTRMPI
# include <mpi.h>
# endif

# include "Qtr.h"
# include "Pointers.h"


using std::cout;
using std::endl;

using namespace QTR_NS;

int main(int narg, char **arg)
{
    if (narg != 2) {
        printf("Syntax: qtr in.file\n");
        exit(1);
    }
    
#ifdef QTRMPI
    MPI_Init(&narg, &arg);
#endif
    
    QTR *qtr;
    
    /* Activate qtr */
    
# ifdef QTRMPI

    qtr = new QTR(MPI_COMM_WORLD,arg);

# else
    qtr = new QTR(arg);
# endif  
    
    qtr->run();
    
    delete qtr;
    
# ifdef  QTRMPI
    MPI_Finalize();
# endif
    
    return 0;  
}


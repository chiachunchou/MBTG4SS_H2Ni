// ==============================================================================
//
//  Constants.cpp
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

#include <string.h>

#include "Constants.h"

const char *elementArray[] =
{
    "Unknown",
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar",
    "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og", NULL};

/* ------------------------------------------------------------------------------- */
    
int symbol2atomicNumber(char const * symbol)
{
    int i = 0;
        
    while (elementArray[i] != NULL) {
            
        if (strcmp(symbol, elementArray[i]) == 0)
        {
            return i;
        }
        i++;
    }
    
    return -1;  // invalid symbol
}
/* ------------------------------------------------------------------------------- */
    
char const * atomicNumber2symbol(int n)
{
    return elementArray[n];
}
/* ------------------------------------------------------------------------------- */

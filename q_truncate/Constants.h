// ==============================================================================
//
//  Constants.h
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

#ifndef QTR_CONSTANTS_H
#define QTR_CONSTANTS_H

/* -------------------------------------------- */
// StdUnit unit system :
// Energy   : eV
// Mass     : Molar Mass
// Time     : StdUnitTime
// Length   : A
// Velocity : A / StdUnitTime
// Force    : eV / A
// Temp     : K
/* -------------------------------------------- */

// Avogadro's Number
#define NA (6.02214078E+23)

// Boltzmann Constant
// unit (J*K^-1)
#define KB   (1.3806488E-23)
// unit (eV*K^-1)
#define KBEV (8.617332384960E-5)

// Gas Constant
#define R_GAS_KMS (8.3144621)

// Bohr Radius
// unit (A)
#define R_BOHR_A (5.291772105126203E-1)
// unit (m)
#define R_BOHR_M (5.291772105126203E-11)

// Electron Charge
// unit (coulomb)
#define EE (1.602176565E-19)

// Electron Mass
// unit (kg)
#define EM_KG (9.10938291E-31)

// Fine Structure Constant
#define FSC (7.29735257E-3)

// Planck's Constant

// unit (J*s)
#define H_PLANCK_JS (6.62606957E-34)
// unit (eV*s)
#define H_PLANCK_EVS (4.135667516E-15)
// unit (eV*fs)
#define H_PLANCK_EVFS (4.135667516)
// unit (eV*StdUnitTime)
#define H_PLANCK_EVSTD (4.06234005040383E-1)

// Plank's hbar

// unit (J*s)
#define HBAR_PLANCK_JS (1.054571725336289E-34)
// unit (eV*s)
#define HBAR_PLANCK_EVS (6.5821192815598E-16)
// unit (eV*fs)
#define HBAR_PLANCK_EVFS (6.5821192815598E-1)
// unit (eV*StdUnitTime)
#define HBAR_PLANCK_EVSTD (6.46541499541949E-2)

// Unit conversion
#define CAL_TO_J (4.18400)
#define WN_TO_HARTREE (4.5563352708321300E-06)

// Speed of Light
// unit (m*s^-1)
#define CC_M (299792458)

// Unit Converter
// unit (StdUnitTime*fs^-1)
#define FS_TO_EV_TIME (0.09822694969282549)

// Length
#define NM_TO_M (1E-9)

// Pi
#define PI (3.1415926535897932384626433832795028841971693993751)

// Functions

int symbol2atomicNumber(char const * symbol);
char const * atomicNumber2symbol(int n);
        

#endif /* QTR_CONSTANTS_H */

#include "pzlog.h"
#include "TPZStiffFracMemory.h"

#include <iostream>
#include <cmath>


TPZStiffFracMemory::TPZStiffFracMemory() : TPZMonoPhasicMemoryDFN() 
{

  fmu = 0.;
  fE = 0.;
  fnu = 0.;
  fSigmaConf = 0.;
  //this->SetCurrentState();
}


TPZStiffFracMemory::~TPZStiffFracMemory()
{
  
}

/** @brief copy constructor */
TPZStiffFracMemory::TPZStiffFracMemory(const TPZStiffFracMemory &copy) : TPZMonoPhasicMemoryDFN(copy) {
    
}

/** @brief operator equal */
TPZStiffFracMemory & TPZStiffFracMemory::operator=(const TPZStiffFracMemory &copy) {
    fmu = copy.fmu;
    fE = copy.fE;
    fnu = copy.fnu;
    fSigmaConf = copy.fSigmaConf;
    
}


void TPZStiffFracMemory::Viscosity(REAL p, REAL &FluidViscosity, REAL &dFluidViscosityDp) const
{
  // Constant for a while.
  FluidViscosity = fmu;
  dFluidViscosityDp = 0;
}



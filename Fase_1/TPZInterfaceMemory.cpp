#include "pzlog.h"
#include "TPZInterfaceMemory.h"

#include <iostream>
#include <cmath>


TPZInterfaceMemory::TPZInterfaceMemory() : TPZMonoPhasicMemoryDFN()
{

  fmu = 0.;
  fE = 0.;
  fnu = 0.;
  fSigmaConf = 0.;
  //this->SetCurrentState();
}


TPZInterfaceMemory::~TPZInterfaceMemory()
{
  
}

/** @brief copy constructor */
TPZInterfaceMemory::TPZInterfaceMemory(const TPZInterfaceMemory &copy) : TPZMonoPhasicMemoryDFN(copy) {
    
}

/** @brief operator equal */
TPZInterfaceMemory & TPZInterfaceMemory::operator=(const TPZInterfaceMemory &copy) {
    fmu = copy.fmu;
    fE = copy.fE;
    fnu = copy.fnu;
    fSigmaConf = copy.fSigmaConf;
    
}


void TPZInterfaceMemory::Viscosity(REAL p, REAL &FluidViscosity, REAL &dFluidViscosityDp) const
{
  // Constant for a while.
  FluidViscosity = fmu;
  dFluidViscosityDp = 0;
}



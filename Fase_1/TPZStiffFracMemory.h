#ifndef TPZStiffFracMemory_h
#define TPZStiffFracMemory_h

#include "pzvec.h"
#include "pzfmatrix.h"
#include "TPZMonoPhasicMemoryDFN.h"

/**
 * @author Omar Duran and Nathan Shauer
 * @since 19/08/2014
 * @brief class to store data of fracture propagation simulation
 */
class TPZStiffFracMemory : public TPZMonoPhasicMemoryDFN {
  
public:
  
    /** @brief Default Constructor */
    TPZStiffFracMemory();
  
    /** @brief Destructor */
    ~TPZStiffFracMemory();
  
    /** @brief copy constructor */
    TPZStiffFracMemory(const TPZStiffFracMemory &copy);
  
    /** @brief operator equal */
    TPZStiffFracMemory &operator=(const TPZStiffFracMemory &copy);
  
private:
  
    /** @brief State: Stiffness or Mass Matrix Calculations */
    enum nState { n = 0, nplusone = 1 };
    nState State;
  
    /** @brief Fluid Viscosity - Pa.s */
    REAL fmu;
  
    /** @brief Fracture height */
    REAL fheight;
  
    /** @brief Reservoir young modulus */
    REAL fE;
  
    /** @brief Reservoir poisson */
    REAL fnu;
  
    /** @brief Confinement Stress */
    REAL fSigmaConf;
    
    /** @brief Simulation time step */
    REAL fDeltaT;
    
    /** @brief Simulation current time */
    REAL fTime;
    
  
public:
    
    /** @brief Enum for materials Ids */
    enum MadIds {EMatDarcy = 1, EBCBottom = 2, EBCRight = 3, EBCTop = 4, EBCLeft = 5, EMatFrac = 6, EMatInterFrac = 20, EBCInlet = 7, EBCOutlet = 8 , EBCAuxBottom = 10};
    
    /** @brief Sets the postprocess file for fracture vtk */
    void SetPostProcessFileName(std::string &postProcessFileName);
  
    /** @brief Returns the postprocess file for fracture vtk */
    std::string PostProcessFileName();
  
    /** @brief Set fluid viscosity. */
    void SetViscosity(REAL mu){this->fmu = mu;}
  
    /** @brief Returns fluid viscosity. */
    REAL Viscosity() const {return this->fmu;}
  
    /**
     * @brief Fluid viscosity function. \f$ FluidViscosity = Visc( p ) \f$
     * @param p pressure
     */
    void Viscosity(REAL p, REAL &FluidViscosity, REAL &dFluidViscosityDp) const;
  
    /** @brief Defines simulation elasticity */
    void SetE(REAL E){ this->fE = E;}
  
    /** @brief Returns simulation elasticity */
    REAL E() const {return this->fE;}
  
    /** @brief Defines simulation poisson */
    void SetPoisson(REAL nu){ this->fnu = nu;}
  
    /** @brief Returns simulation poisson */
    REAL Poisson() const {return this->fnu;}
  
    /** @brief Defines simulation Confinement stress */
    void SetSigmaConf(REAL SigmaConf){ this->fSigmaConf = SigmaConf;}
  
    /** @brief Returns simulation confinement stress */
    REAL SigmaConf() const {return this->fSigmaConf;}
    
    /** @brief Defines simulation time step. */
    void SetTimeStep(REAL timestep){ this->fDeltaT = timestep;}
    
    /** @brief Returns simulation time step. */
    REAL TimeStep() const {return this->fDeltaT;}

    
};

#endif

//
//  TPZMatFractureBB.h
//  Benchmark0a
//
//  Created by Pablo Carvalho on 20/08/18.
//

#ifndef TPZMatFractureBB_h
#define TPZMatFractureBB_h

#include <stdio.h>
#include <iostream>
#include "pzdiscgal.h"
#include "TPZMatWithMem.h"
#include "tpzautopointer.h"
#include "TPZMemoryFracDFN.h"
#include "TPZSimulationData.h"

/// Material which implements a Lagrange Multiplier

template <class TMEM = TPZMemoryFracDFN>
class TPZMatFractureBB : public TPZMatWithMem<TMEM>
{
    protected:
    
    // Number of state variables
    int fNStateVariables;
    
    // Dimensiona associated with the material
    int fDimension;
    
    STATE fMultiplier;
    
    // Dimensiona associated with the material
    STATE fFracHsize;
    
    /// Pointer of Simulation data
    TPZSimulationData * m_simulation_data;
    
    
    public :
    // Simple constructor
    TPZMatFractureBB() : TPZRegisterClassId(&TPZMatFractureBB::ClassId), TPZMatWithMem<TMEM>()
    {
        
    }
    // Constructor with the index of the material object within the vector
    TPZMatFractureBB(int nummat, int dimension, int nstate) : TPZRegisterClassId(&TPZMatFractureBB::ClassId),TPZMatWithMem<TMEM>(nummat), fNStateVariables(nstate), fDimension(dimension), fMultiplier(1.)
    {
        m_simulation_data = NULL;
    }
    
    // Copy constructor
    TPZMatFractureBB(const TPZMatFractureBB &copy) : TPZRegisterClassId(&TPZMatFractureBB::ClassId),TPZMatWithMem<TMEM>(copy), fNStateVariables(copy.fNStateVariables), fDimension(copy.fDimension), fMultiplier(copy.fMultiplier)
    {
        m_simulation_data = copy.m_simulation_data;
    }
    
    TPZMatFractureBB &operator=(const TPZMatFractureBB &copy)
    {
        if(&copy == this){
            return *this;
        }
        fNStateVariables = copy.fNStateVariables;
        fDimension = copy.fDimension;
        fMultiplier = copy.fMultiplier;
        return *this;
    }
    
    
    // Destructor
    ~TPZMatFractureBB()
    {
        
    }
    
    // Returns the integrable dimension of the material
    int Dimension() const
    {
        return fDimension;
    }
    
    void SetMultiplier(STATE mult)
    {
        fMultiplier = mult;
    }

    void SetFracHSize(STATE fracHsize)
    {
        fFracHsize = fracHsize;
    }

    
    STATE GetFracHSize()
    {
        return fFracHsize;
    }
    
    virtual std::string Name()
    {
        return "TPZMatFractureBB";
    }
    
    // Fill material data parameter with necessary requirements for the ContributeInterface method.
    // Here, in base class, all requirements are considered as necessary. \n
    // Each derived class may optimize performance by selecting only the necessary data.
    void FillDataRequirementsInterface(TPZMaterialData &data)
    {
        data.SetAllRequirements(false);
    }
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef);
    
    void ContributeBC(TPZMaterialData &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    int NStateVariables()
    {
        return fNStateVariables;
    }
    
    // Updates the leak off memory
    void UpdateMemory(TPZVec<TPZMaterialData> &datavec);
    
    // Updates the leak off memory
    void UpdateMemory(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavec);
    

public:

    // Saves the element data to a stream
    void Write(TPZStream &buf, int withclassid);
    
    
    // Reads the element data from a stream
    void Read(TPZStream &buf, void *context);
    

};


#endif /* TPZMatFractureBB_h */

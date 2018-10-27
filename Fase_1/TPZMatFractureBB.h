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

/// Material which implements the Barton-Bandis fracture opnenig model

template <class TMEM = TPZMemoryFracDFN>
class TPZMatFractureBB : public TPZMatWithMem<TMEM>
{
    protected:
    
    // Number of state variables
    int fNStateVariables;
    
    // Dimensiona associated with the material
    int fDimension;
    
    STATE fMultiplier;
    
    // Initial fracture opening
    STATE m_a0;
    
    // Maximum fracture opening
    STATE m_Vm;
    
    // Initial normal stiffness
    STATE m_Kni;
    
    
    /// Pointer of Simulation data
    TPZSimulationData * m_simulation_data;
    
    
    public :
    // Simple constructor
    TPZMatFractureBB() : TPZRegisterClassId(&TPZMatFractureBB::ClassId), TPZMatWithMem<TMEM>(), fNStateVariables(0.), fDimension(0.), fMultiplier(0.), m_a0(0.), m_Vm(0.), m_Kni(0.)
    {
        m_simulation_data = NULL;
    }
    // Constructor with the index of the material object within the vector
    TPZMatFractureBB(int nummat, int dimension, int nstate) : TPZRegisterClassId(&TPZMatFractureBB::ClassId),TPZMatWithMem<TMEM>(nummat), fNStateVariables(nstate), fDimension(dimension), fMultiplier(1.), m_a0(0.), m_Vm(0.), m_Kni(0.)
    {
        m_simulation_data = NULL;
    }
    
    // Copy constructor
    TPZMatFractureBB(const TPZMatFractureBB &copy) : TPZRegisterClassId(&TPZMatFractureBB::ClassId),TPZMatWithMem<TMEM>(copy), fNStateVariables(copy.fNStateVariables), fDimension(copy.fDimension), fMultiplier(copy.fMultiplier), m_a0(copy.m_a0), m_Vm(copy.m_Vm), m_Kni(copy.m_Kni)
    {
        m_simulation_data = copy.m_simulation_data;
    }
    
    /// Set the pointer of Simulation data object
    void SetSimulationData(TPZSimulationData * simulation_data){
        m_simulation_data = simulation_data;
    }
    
    TPZMatFractureBB &operator=(const TPZMatFractureBB &copy)
    {
        if(&copy == this){
            return *this;
        }
        fNStateVariables = copy.fNStateVariables;
        fDimension = copy.fDimension;
        fMultiplier = copy.fMultiplier;
        m_a0 = copy.m_a0;
        m_Vm = copy.m_Vm;
        m_Kni = copy.m_Kni;
        
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
    
    void Set_a0(STATE a0)
    {
        m_a0 = a0;
    }

    STATE Get_a0()
    {
        return m_a0;
    }

    void Set_Vm(STATE Vm)
    {
        m_Vm = Vm;
    }
    
    STATE Get_Vm()
    {
        return m_Vm;
    }
    
    void Set_Kni(STATE Kni)
    {
        m_Kni = Kni;
    }
    
    STATE Get_Kni()
    {
        return m_Kni;
    }
    
    virtual std::string Name()
    {
        return "TPZMatFractureBB";
    }
    
    // Fill material data parameter with necessary requirements
    void FillDataRequirements(TPZMaterialData &data)
    {
        data.SetAllRequirements(false);
        data.fNeedsNormal = true;
        data.fNeedsNeighborSol = true;
        data.fNeedsSol = true;
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
    
    void SetNStateVariables(int nstateVariables)
    {
        fNStateVariables = nstateVariables;
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

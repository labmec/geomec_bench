//
//  TPZLagrangeInterface.h
//  Benchmark0a
//
//  Created by Pablo Carvalho on 20/08/18.
//

#ifndef TPZLagrangeInterface_h
#define TPZLagrangeInterface_h

#include <stdio.h>
#include <iostream>
#include "pzdiscgal.h"
#include "TPZMatWithMem.h"
#include "tpzautopointer.h"
#include "TPZInterfaceMemory.h"
#include "TPZSimulationData.h"

/// Material which implements a Lagrange Multiplier

template <class TMEM = TPZInterfaceMemory>
class TPZLagrangeInterface : public TPZMatWithMem< TMEM, TPZDiscontinuousGalerkin>
{
    protected:
    
    // Number of state variables
    int fNStateVariables;
    
    // Dimensiona associated with the material
    int fDimension;
    
    STATE fMultiplier;
    
    /// Pointer of Simulation data
    TPZSimulationData * m_simulation_data;
    
    public :
    // Simple constructor
    TPZLagrangeInterface() : TPZRegisterClassId(&TPZLagrangeInterface::ClassId), TPZMatWithMem<TMEM,
    TPZDiscontinuousGalerkin >()
    {
        
    }
    // Constructor with the index of the material object within the vector
    TPZLagrangeInterface(int nummat, int dimension, int nstate) : TPZRegisterClassId(&TPZLagrangeInterface::ClassId),TPZMatWithMem<TMEM,
    TPZDiscontinuousGalerkin>(nummat), fNStateVariables(nstate), fDimension(dimension), fMultiplier(1.)
    {
        
    }
    
    // Copy constructor
    TPZLagrangeInterface(const TPZLagrangeInterface &copy) : TPZRegisterClassId(&TPZLagrangeInterface::ClassId),TPZMatWithMem<TMEM,
    TPZDiscontinuousGalerkin>(copy), fNStateVariables(copy.fNStateVariables), fDimension(copy.fDimension), fMultiplier(copy.fMultiplier)
    {
        
    }
    
    TPZLagrangeInterface &operator=(const TPZLagrangeInterface &copy)
    {
        TPZDiscontinuousGalerkin::operator=(copy);
        fNStateVariables = copy.fNStateVariables;
        fDimension = copy.fDimension;
        fMultiplier = copy.fMultiplier;
        return *this;
    }
    
//    TPZLagrangeInterface *NewMaterial()
//    {
//        return new TPZLagrangeInterface(*this);
//    }
    
    // Destructor
    ~TPZLagrangeInterface()
    {
        
    }
    
    // Returns the integrable dimension of the material
    virtual int Dimension() const
    {
        return fDimension;
    }
    
    virtual void SetMultiplier(STATE mult)
    {
        fMultiplier = mult;
    }
    
    
    virtual std::string Name()
    {
        return "TPZLagrangeInterface";
    }
    
    /// Set the pointer of Simulation data object
    void SetSimulationData(TPZSimulationData * simulation_data){
        m_simulation_data = simulation_data;
    }
    
    // Fill material data parameter with necessary requirements for the ContributeInterface method.
    // Here, in base class, all requirements are considered as necessary. \n
    // Each derived class may optimize performance by selecting only the necessary data.
    virtual void FillDataRequirementsInterface(TPZMaterialData &data)
    {
        data.SetAllRequirements(false);
        data.fNeedsNeighborSol = true;
    }

    virtual void FillDataRequirements(TPZMaterialData &data)
    {
        data.SetAllRequirements(false);
        data.fNeedsNeighborSol = true;
    }
    
    // Contribute methods
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    // It computes a contribution to stiffness matrix and load vector at one integration point
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    // Computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    // It computes a contribution to residual vector at one integration point
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef);
    
    // Computes a contribution to residual vector at one integration point
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ef)
    {
        DebugStop();
        ContributeInterface(data, dataleft[0], dataright[0], weight, ef);
    }
    
    virtual int NStateVariables() const
    {
        return fNStateVariables;
    }
    
    // Updates the leak off memory
    void UpdateMemory(TPZVec<TPZMaterialData> &datavec);
    
    // Updates the leak off memory
    void UpdateMemory(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavec);
    
    // Inner vec
    STATE InnerVec(TPZFMatrix<STATE>  &S, TPZManVector<STATE,3>  &T);

    STATE InnerVec(TPZManVector<STATE,3>  &S, TPZManVector<STATE,3>  &T);
    
public:
    
    // Unique identifier for serialization purposes
    virtual int ClassId() const;
    
    
    // Saves the element data to a stream
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    
    // Reads the element data from a stream
    virtual void Read(TPZStream &buf, void *context);
    

};


#endif /* TPZLagrangeInterface_h */

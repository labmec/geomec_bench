/*
 *  TPZDarcy2DMaterialMem.h
 *  PZ
 *
 *  Created by Pablo G S Carvalho on 08/09/16.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef TPZDarcy2DMaterialMem_h
#define TPZDarcy2DMaterialMem_h

#include <stdio.h>
#include "TPZMatWithMem.h"
#include "tpzautopointer.h"
#include "pzdiscgal.h"
#include "pzfmatrix.h"
#include "pzbndcond.h"
#include "pzlog.h"
#include "tpzautopointer.h"
#include "TPZMaterial.h"
#include "TPZElastoPlasticMemoryDFN.h"
#include "pzporoelastoplasticmem.h"
#include "TPZMonoPhasicMemoryDFN.h"
#include "TPZMemoryBCDFN.h"
#include "TPZMemoryDFN.h"
#include "TPZMemoryFracDFN.h"
#include "TPZSimulationData.h"


template <class TMEM>
class TPZDarcy2DMaterialMem : public TPZMatWithMem<TMEM>  {
    
//protected:
    
    /// dimension of the material
    int fDimension;
    
    /// Aproximation Space for velocity
    int fSpace;
    
    /** @brief Forcing function value */
    REAL ff;
    
    /// viscosidade
    STATE fViscosity;
    
    /** @brief Medium permeability. Coeficient which multiplies the gradient operator*/
    REAL fk;

    /** @brief Medium vertical permeability. Coeficient which multiplies the gradient operator*/
    REAL fkv;
    
    /** @brief Medium hotizontal permeability. Coeficient which multiplies the gradient operator*/
    REAL fkh;
    
    
    /** @brief Medium porosity.*/
    REAL fPhi;
    
    /** @brief permeability tensor. Coeficient which multiplies the gradient operator*/
    TPZFNMatrix<9,REAL> fTensorK;
    
    /** @brief inverse of the permeability tensor.*/
    TPZFNMatrix<9,REAL> fInvK;
    
    /// termo contrario a beta na sua formulacao (para ser conforme a literatura)
    STATE fTheta;
    
    /** @brief Pointer to forcing function, it is the Permeability and its inverse */
    TPZAutoPointer<TPZFunction<STATE> > fPermeabilityFunction;
    
    /// Pointer of Simulation data
    TPZSimulationData * m_simulation_data;
    
public:
    
    
    REAL BIGNUMBER = TPZMaterial::gBigNumber;
    
    /**
     * Empty Constructor
     */
    TPZDarcy2DMaterialMem();
    
    /** Creates a material object and inserts it in the vector of
     *  material pointers of the mesh.
     */
    TPZDarcy2DMaterialMem(int matid, int dimension, int space, STATE theta);
    
    
    /** Creates a material object based on the referred object and
     *  inserts it in the vector of material pointers of the mesh.
     */
    TPZDarcy2DMaterialMem(const TPZDarcy2DMaterialMem &mat);
    
    TPZDarcy2DMaterialMem &operator=(const TPZDarcy2DMaterialMem & other);
    
    /**
     * Destructor
     */
    ~TPZDarcy2DMaterialMem();
    
    /** Fill material data parameter with necessary requirements for the
     * Contribute method. Here, in base class, all requirements are considered
     * as necessary. Each derived class may optimize performance by selecting
     * only the necessary data.
     * @since April 10, 2007
     */
    
    void FillDataRequirements(TPZMaterialData &data);
    
    void FillDataRequirements(TPZVec<TPZMaterialData> &datavec);
    
    virtual void FillDataRequirementsInterface(TPZMaterialData &data);
    
    void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec);
    
    void SetPermeabilityXY(TPZFMatrix<REAL> perm) {
        fkh = perm(0,0);
        fkv = perm(1,1);
        fTensorK.Zero();
        fInvK.Zero();
        for (int i=0; i<fDimension; i++) {
            fTensorK(i,i) = perm(i,i);
            fInvK(i,i) = 1./perm(i,i);
        }
    }
    
    void SetPermeability(REAL perm) {
        fk = perm;
        fTensorK.Zero();
        fInvK.Zero();
        for (int i=0; i<fDimension; i++) {
            fTensorK(i,i) = perm;
            fInvK(i,i) = 1./perm;
        }
    }

    REAL GetPermeability(){
        return fk;
    }
    
    //Set the permeability tensor and inverser tensor
    void SetPermeabilityTensor(TPZFMatrix<REAL> K, TPZFMatrix<REAL> invK){
        
        //       if(K.Rows() != fDim || K.Cols() != fDim) DebugStop();
        //       if(K.Rows()!=invK.Rows() || K.Cols()!=invK.Cols()) DebugStop();
        
        fTensorK = K;
        fInvK = invK;
    }
    
    void SetViscosity(REAL visc) {
        fViscosity = visc;
    }

    void SetPorosity(REAL poro) {
        fPhi = poro;
    }
    
    REAL GetPorosity() {
        return fPhi;
    }
    
    void SetInternalFlux(REAL flux) {
        ff = flux;
    }
    
    /// Set the pointer of Simulation data object
    void SetSimulationData(TPZSimulationData * simulation_data){
        m_simulation_data = simulation_data;
    }
    
    /** returns the name of the material */
    std::string Name() {
        return "TPZDarcy2DMaterial";
    }
    
    /** returns the integrable dimension of the material */
    int Dimension() const {return fDimension;}
    
    /** returns the number of state variables associated with the material */
    virtual int NStateVariables() const {return 1;} // for hdiv are 3, plus pressure, so 3 + 1 = 4 itapopo
    
    /** print out the data associated with the material */
    void Print(std::ostream &out = std::cout);
    
    /** returns the variable index associated with the name */
    int VariableIndex(const std::string &name);
    
    /** returns the number of variables associated with the variable
     indexed by var.  var is obtained by calling VariableIndex */
    int NSolutionVariables(int var);
    
    /** Computes the divergence over the parametric space */
    void ComputeDivergenceOnDeformed(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi);
    
    void ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi);
    
    
    void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
    
    /** returns the solution associated with the var index based on
     * the finite element approximation */
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    
    /** index of velocity */
    int VIndex(){ return 0; }
    
    /** index of pressure */
    int PIndex(){ return 1; }
    
    /** inner product of two tensors. See Gurtin (2003), p. 5. */
    template <typename TVar>
    TVar Inner(TPZFMatrix<TVar> &S, TPZFMatrix<TVar> &T);
    
    /** inner product of two vectors. See Gurtin (2003), p. 5. */
    template <typename TVar>
    TVar InnerVec(TPZFMatrix<TVar> &S, TPZFMatrix<TVar> &T);
    
    /** trace of the tensor GradU = Div(U)*/
    STATE Tr(TPZFMatrix<REAL> &GradU );
    
    /** transpose of the tensor GradU = Div(U)*/
    STATE Transpose(TPZFMatrix<REAL> &GradU );
    
    /** Fill the vector of gradient for each phi */
    void FillGradPhi(TPZMaterialData &dataV, TPZVec< TPZFMatrix<REAL> > &GradPhi);
    
    /// transform a H1 data structure to a vector data structure
    void FillVecShapeIndex(TPZMaterialData &data);
    
    /** @brief Not used contribute methods */
    virtual void Contribute(TPZMaterialData &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){DebugStop();}
    
    // Contribute Methods being used - Multiphysics
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
    {
        DebugStop();
    }
    
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    
    // Contribute Methods being used
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    /**
     * Save the element data to a stream
     */
    void Write(TPZStream &buf, int withclassid);
    
    /**
     * Read the element data from a stream
     */
    void Read(TPZStream &buf, void *context);
    
    
    virtual int NEvalErrors() {return 6;}
    
    /**
     * It computes errors contribution in differents spaces.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     */
    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors);
    
    
    void porosity(long gp_index, REAL &phi_n);
    
};

#endif

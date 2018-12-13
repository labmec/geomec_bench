//
//  TPZMemoryDFN.cpp
//  Benchmark0a
//
//  Created by Pablo Carvalho on 17/09/18.
//


#ifndef TPZMemoryDFN_h
#define TPZMemoryDFN_h

#include <stdio.h>
#include "TPZElastoPlasticMemoryDFN.h"
#include "TPZMonoPhasicMemoryDFN.h"

class TPZMemoryDFN : public TPZMonoPhasicMemoryDFN, public TPZElastoPlasticMemoryDFN {
    
private:
    
    /// Biot-Willis coefficient
    REAL m_alpha;
    
    // Point coordenates
    TPZManVector<STATE,3> m_coord;
    
    
    
public:
    
    /// Default constructor
    TPZMemoryDFN();
    
    /// Copy constructor
    TPZMemoryDFN(const TPZMemoryDFN & other);
    
    /// Assignement constructor
    const TPZMemoryDFN & operator=(const TPZMemoryDFN & other);
    
    /// Desconstructor
    virtual ~TPZMemoryDFN();
    
    /// Class name
    const std::string Name() const;
    
    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// Read class attributes
    virtual void Read(TPZStream &buf, void *context);
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    /// Print class attributes
    friend std::ostream & operator<<( std::ostream& out, const TPZMemoryDFN & memory ){
        memory.Print(out);
        return out;
    }
    
    REAL Permeability(REAL k0, REAL phi0, REAL nu, REAL E){
        
        REAL Sigma_kk_ef = GetSigma_n()[0]+GetSigma_n()[3]+GetSigma_n()[5];
        
        REAL strain_v_n = (1.-2.*nu)*(Sigma_kk_ef)/E + 0.0013171510898339536;
        
        //REAL strain_v_n = (1.-2.*nu)*(Sigma_kk_ef)/E+0.0013171562035842917;
  
        //strain_v_n = 0.00000798717690920986;
        REAL phi   = 1. - (1. - phi0) * exp(-strain_v_n);
        //REAL phi   = 1. - (1. - phi0) * exp(u_v);
        REAL varphi = phi/phi0;
 
        REAL perm = k0*pow(varphi, 60.);
        
        return perm;
    }
    
    virtual int ClassId() const;
    
    void SetAlpha(REAL alpha){
        m_alpha  = alpha;
    }
    
    REAL GetAlpha(){
        return m_alpha;
    }
    
    /// Set the coordinates
    void SetCoord(TPZManVector<STATE,3> & Xcoord){
        m_coord = Xcoord;
    }
    
    /// Get the coordinates
    TPZManVector<STATE,3> GetCoord(){
        return m_coord;
    }
    
};

#endif /* TPZMemoryDFN_h */


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
    
    REAL Permeability(STATE k0, STATE phi0, STATE nu, STATE E){
        
        REAL sigma_v_n = GetSigma_n().I1()/3.;
        REAL strain_v_n = (1.-2.*nu)*sigma_v_n/E;
        REAL phi   = 1. - (1. - phi0) * exp(strain_v_n);
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
    
};

#endif /* TPZMemoryDFN_h */


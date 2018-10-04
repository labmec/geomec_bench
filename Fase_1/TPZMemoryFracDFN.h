//
//  TPZMemoryFracDFN.cpp
//  Benchmark0a
//
//  Created by Pablo Carvalho on 17/09/18.
//


#ifndef TPZMemoryFracDFN_h
#define TPZMemoryFracDFN_h

#include <stdio.h>
#include "TPZElastoPlasticMemoryDFN.h"
#include "TPZMonoPhasicMemoryDFN.h"

class TPZMemoryFracDFN : public TPZMonoPhasicMemoryDFN, public TPZElastoPlasticMemoryDFN {
    
private:
    
    /// Biot-Willis coefficient
    REAL m_alpha;
    
    /// Relative displacement
    TPZManVector<REAL,3> m_uR;
    
public:
    
    /// Default constructor
    TPZMemoryFracDFN();
    
    /// Copy constructor
    TPZMemoryFracDFN(const TPZMemoryFracDFN & other);
    
    /// Assignement constructor
    const TPZMemoryFracDFN & operator=(const TPZMemoryFracDFN & other);
    
    /// Desconstructor
    virtual ~TPZMemoryFracDFN();
    
    /// Class name
    const std::string Name() const;
    
    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// Read class attributes
    virtual void Read(TPZStream &buf, void *context);
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    /// Print class attributes
    friend std::ostream & operator<<( std::ostream& out, const TPZMemoryFracDFN & memory ){
        memory.Print(out);
        return out;
    }
    
    virtual int ClassId() const;
    
    void SetAlpha(REAL alpha){
        m_alpha  = alpha;
    }
    
    REAL GetAlpha(){
        return m_alpha;
    }
    
    /// Set the relative displacement
    void SetuR_n(TPZManVector<REAL,3> & uR){
        m_uR = uR;
    }
    
    /// Get the relative displacement
    TPZManVector<REAL,3> & GetuR(){
        return m_uR;
    }

    
};

#endif /* TPZMemoryFracDFN_h */


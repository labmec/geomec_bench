//
//  TPZMemoryBCDFN.cpp
//  Benchmark0a
//
//  Created by Pablo Carvalho on 17/09/18.
//


#ifndef TPZMemoryBCDFN_h
#define TPZMemoryBCDFN_h

#include <stdio.h>
#include "TPZElastoPlasticMemoryBCDFN.h"
#include "TPZMonoPhasicMemoryBCDFN.h"

class TPZMemoryBCDFN : public TPZMonoPhasicMemoryBCDFN, public TPZElastoPlasticMemoryBCDFN {
    
private:
    
    /// Biot-Willis coefficient
    REAL m_alpha;
    
    
public:
    
    /// Default constructor
    TPZMemoryBCDFN();
    
    /// Copy constructor
    TPZMemoryBCDFN(const TPZMemoryBCDFN & other);
    
    /// Assignement constructor
    const TPZMemoryBCDFN & operator=(const TPZMemoryBCDFN & other);
    
    /// Desconstructor
    virtual ~TPZMemoryBCDFN();
    
    /// Class name
    const std::string Name() const;
    
    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// Read class attributes
    virtual void Read(TPZStream &buf, void *context);
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    /// Print class attributes
    friend std::ostream & operator<<( std::ostream& out, const TPZMemoryBCDFN & memory ){
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
    

    
};

#endif /* TPZMemoryBCDFN_h */


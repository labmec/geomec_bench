//
//  TPZMemoryFracDFN.cpp
//  Benchmark0a
//
//  Created by Pablo Carvalho on 17/09/18.
//


#ifndef TPZMemoryFracDFN_h
#define TPZMemoryFracDFN_h

#include <stdio.h>
#include "TPZElastoPlasticMemoryFracDFN.h"
#include "TPZMonoPhasicMemoryFracDFN.h"

class TPZMemoryFracDFN : public TPZMonoPhasicMemoryFracDFN, public TPZElastoPlasticMemoryFracDFN {
    
private:
    
    /// Biot-Willis coefficient
    REAL m_alpha;
    
    /// Relative displacement (fracture oppening)
    TPZManVector<REAL,3> m_uR;
    
    // Point coordenates
    TPZManVector<REAL,3> m_coord;
    
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

    /// Set the coordinates
    void SetCoord(TPZManVector<REAL,3> & Xcoord){
        m_coord = Xcoord;
    }
    
    /// Get the coordinates
    TPZManVector<REAL,3> GetCoord(){
        return m_coord;
    }
    
};

#endif /* TPZMemoryFracDFN_h */


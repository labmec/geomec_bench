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
    STATE m_alpha;

    /// Initial fracture closure
    STATE m_Du_0;

    /// Last fracture closure
    TPZManVector<STATE,3> m_Du;
    
    /// Current fracture closure
    TPZManVector<STATE,3> m_Du_n;
    
    // Point coordenates
    TPZManVector<STATE,3> m_coord;
    
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
    
    void SetAlpha(STATE alpha){
        m_alpha  = alpha;
    }
    
    STATE GetAlpha(){
        return m_alpha;
    }
    
    /// Set the initial fracture closure
    void SetDu_0(STATE uR){
        m_Du_0 = uR;
    }
    
    /// Get the initial fracture closure
    STATE & GetDu_0(){
        return m_Du_0;
    }
    
    /// Set the last fracture closure
    void SetDu(TPZManVector<STATE,3> & uR){
        m_Du = uR;
    }
    
    /// Get the last fracture closure
    TPZManVector<STATE,3> & GetDu(){
        return m_Du;
    }
    
    /// Set the current fracture closure
    void SetDu_n(TPZManVector<STATE,3> & uR){
        m_Du_n = uR;
    }
    
    /// Get the current fracture closure
    TPZManVector<STATE,3> & GetDu_n(){
        return m_Du_n;
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

#endif /* TPZMemoryFracDFN_h */


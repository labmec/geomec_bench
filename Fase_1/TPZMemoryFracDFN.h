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
#include "TPZElastoPlasticMemoryDFN.h"

class TPZMemoryFracDFN : public TPZMonoPhasicMemoryFracDFN, public TPZElastoPlasticMemoryFracDFN {
    
private:
    
    /// Biot-Willis coefficient
    STATE m_alpha;

    /// Maximum opening of fracture
    STATE m_Vm;
    
    /// Initial fracture closure
    STATE m_Du_0;

    /// Last fracture closure
    STATE m_Du;
    
    /// Current fracture closure
    STATE m_Du_n;
    
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
    
    REAL Permeability(STATE k0, STATE phi0, STATE nu, STATE E){
        
        if(m_Du_0==0){
            return k0;
        }else{
            
            REAL h = m_Vm - m_Du_n;
            REAL h_0 = m_Vm - m_Du_0;
            //return k0*exp(-268.*(h - h_0)/h_0);
            return k0;
        }

    }
    
    REAL Permeability(REAL k0){
        
        if(m_Du_0==0){
            return k0;
        }else{
 
            
            REAL h = m_Vm - m_Du_n;
            REAL h_0 = m_Vm - m_Du_0;
            REAL perm = k0*exp(268.*(h-h_0)/h_0);
            return k0;
            //return k0;
        }
        
    }
    
    
    virtual int ClassId() const;
    
    void SetVm(STATE Vm){
        m_Vm  = Vm;
    }
    
    STATE GetVm(){
        return m_Vm;
    }
    
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
    void SetDu(STATE uR){
        m_Du = uR;
    }
    
    /// Get the last fracture closure
    STATE & GetDu(){
        return m_Du;
    }
    
    /// Set the current fracture closure
    void SetDu_n(STATE & uR){
        m_Du_n = uR;
    }
    
    /// Get the current fracture closure
    STATE & GetDu_n(){
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


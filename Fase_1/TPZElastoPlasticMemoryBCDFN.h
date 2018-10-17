//
//  TPZElastoPlasticMemoryBCDFN.h
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#ifndef TPZElastoPlasticMemoryBCDFN_h
#define TPZElastoPlasticMemoryBCDFN_h

#include <stdio.h>
#include <iostream>
#include <string>
#include "TPZSavable.h"
#include "Hash/TPZHash.h"
#include "TPZStream.h"
#include "TPZPersistenceManager.h"
#include "TPZTensor.h"
#include "TPZPlasticState.h"
#include "TPZBndCondWithMem.h"

class TPZElastoPlasticMemoryBCDFN{
    
private:
    
    /// Last displacement field
    TPZManVector<REAL,3> m_u;

    /// Current displacement field
    TPZManVector<REAL,3> m_u_n;
    
    /// Last stress state
    TPZTensor<REAL> m_sigma;
    
    /// Current stress state
    TPZTensor<REAL> m_sigma_n;
    
public:
    
    /// Default constructor
    TPZElastoPlasticMemoryBCDFN();
    
    /// Copy constructor
    TPZElastoPlasticMemoryBCDFN(const TPZElastoPlasticMemoryBCDFN & other);
    
    /// Assignement constructor
    const TPZElastoPlasticMemoryBCDFN & operator=(const TPZElastoPlasticMemoryBCDFN & other);
    
    /// Desconstructor
    virtual ~TPZElastoPlasticMemoryBCDFN();
    
    /// Class name
    const std::string Name() const;
    
    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// Read class attributes
    virtual void Read(TPZStream &buf, void *context);
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    /// Print class attributes
    friend std::ostream & operator<<( std::ostream& out, const TPZElastoPlasticMemoryBCDFN & memory ){
        memory.Print(out);
        return out;
    }
    
    virtual int ClassId() const;

    /// Set the last displacement field
    void Setu(TPZManVector<REAL,3> & u){
        m_u = u;
    }
    
    /// Get the last displacement field
    TPZManVector<REAL,3> & Getu(){
        return m_u;
    }

    /// Set the current displacement field
    void Setu_n(TPZManVector<REAL,3> & u_n){
        m_u_n = u_n;
    }
    
    /// Get the current displacement field
    TPZManVector<REAL,3> & Getu_n(){
        return m_u_n;
    }
    
    /// Set the last stress state
    void SetSigma(TPZTensor<REAL> & sigma){
        m_sigma = sigma;
    }
    
    /// Get the last stress state
    TPZTensor<REAL> & GetSigma(){
        return m_sigma;
    }
    
    /// Set the current stress state
    void SetSigma_n(TPZTensor<REAL> & sigma_n){
        m_sigma_n = sigma_n;
    }
    
    /// Get the current stress state
    TPZTensor<REAL> & GetSigma_n(){
        return m_sigma_n;
    }
    
    
};

#endif /* TPZElastoPlasticMemoryBCDFN_h */

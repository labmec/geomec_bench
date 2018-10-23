//
//  TPZElastoPlasticMemoryFracDFN.h
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#ifndef TPZElastoPlasticMemoryFracDFN_h
#define TPZElastoPlasticMemoryFracDFN_h

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

class TPZElastoPlasticMemoryFracDFN{
    
private:
    
    /// Current normal stress state in fractures
    TPZManVector<REAL,3> m_forceFrac_n;

    /// Initial normal stress state in fractures
    TPZManVector<REAL,3> m_forceFrac_0;
    
    /// Current displacement field
    TPZManVector<REAL,3> m_u_n;
    
    /// Last displacement field
    TPZManVector<REAL,3> m_u;
    
    
public:
    
    /// Default constructor
    TPZElastoPlasticMemoryFracDFN();
    
    /// Copy constructor
    TPZElastoPlasticMemoryFracDFN(const TPZElastoPlasticMemoryFracDFN & other);
    
    /// Assignement constructor
    const TPZElastoPlasticMemoryFracDFN & operator=(const TPZElastoPlasticMemoryFracDFN & other);
    
    /// Desconstructor
    virtual ~TPZElastoPlasticMemoryFracDFN();
    
    /// Class name
    const std::string Name() const;
    
    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// Read class attributes
    virtual void Read(TPZStream &buf, void *context);
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    /// Print class attributes
    friend std::ostream & operator<<( std::ostream& out, const TPZElastoPlasticMemoryFracDFN & memory ){
        memory.Print(out);
        return out;
    }
    
    virtual int ClassId() const;
    
    
    /// Set the current stress state
    void SetForceFrac_n(TPZManVector<REAL,3> & ForceFrac_n){
        m_forceFrac_n = ForceFrac_n;
    }
    
    /// Get the current stress state
    TPZManVector<REAL,3> & GetForceFrac_n(){
        return m_forceFrac_n;
    }
    
    /// Set the current displacement field
    void Setu_n(TPZManVector<REAL,3> & u_n){
        m_u_n = u_n;
    }
    
    /// Get the current displacement field
    TPZManVector<REAL,3> & Getu_n(){
        return m_u_n;
    }
    
    /// Set the last displacement field
    void Setu(TPZManVector<REAL,3> & u){
        m_u = u;
    }
    
    /// Get the last displacement field
    TPZManVector<REAL,3> & Getu(){
        return m_u;
    }
    
};

#endif /* TPZElastoPlasticMemoryFracDFN_h */

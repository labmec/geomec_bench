//
//  TPZMonoPhasicMemoryDFN.h
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#ifndef TPZMonoPhasicMemoryDFN_h
#define TPZMonoPhasicMemoryDFN_h

#include <stdio.h>
#include <iostream>
#include <string>
#include "TPZSavable.h"
#include "Hash/TPZHash.h"
#include "TPZStream.h"
#include "TPZPersistenceManager.h"
#include "pzfmatrix.h"

class TPZMonoPhasicMemoryDFN{
    
private:
    
    /// pore pressure at initial REAL
    REAL m_p_0;
    
    /// pore pressure
    REAL m_p;
    
    /// pore pressure
    REAL m_p_n;
    
    /// absolute permeability at initial REAL
    TPZFNMatrix<9,REAL> m_kappa_0;
    
    /// absolute permeability
    TPZFNMatrix<9,REAL> m_kappa;
    
    /// lagrangian porosity at intial REAL
    REAL m_phi_0;
    
    /// lagrangian porosity
    REAL m_phi;
    
public:
    
    /// Default constructor
    TPZMonoPhasicMemoryDFN();
    
    /// Copy constructor
    TPZMonoPhasicMemoryDFN(const TPZMonoPhasicMemoryDFN & other);
    
    /// Assignement constructor
    const TPZMonoPhasicMemoryDFN & operator=(const TPZMonoPhasicMemoryDFN & other);
    
    /// Desconstructor
    virtual ~TPZMonoPhasicMemoryDFN();
    
    /// Class name
    const std::string Name() const;
    
    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// Read class attributes
    virtual void Read(TPZStream &buf, void *context);
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    /// Print class attributes
    friend std::ostream & operator<<( std::ostream& out, const TPZMonoPhasicMemoryDFN & memory ){
        memory.Print(out);
        return out;
    }
    
    virtual int ClassId() const;
    
    
    /// set and get methods
    
    /// Set pore pressure at initial REAL
    void Setp_0(REAL p_0)
    {
        m_p_0 = p_0;
    }
    
    /// Get pore pressure at initial REAL
    REAL p_0()
    {
        return m_p_0;
    }
    
    
    /// Set pore pressure at current REAL
    void Setp(REAL p)
    {
        m_p = p;
    }
    
    /// Get pore pressure at current REAL
    REAL p()
    {
        return m_p;
    }
    
    
    /// Set pore pressure at last REAL
    void Setp_n(REAL p_n)
    {
        m_p_n = p_n;
    }
    
    /// Get pore pressure at last REAL
    virtual REAL p_n()
    {
        return m_p_n;
    }
    
    
    /// Set absolute permeability at initial REAL
    void Setkappa_0(TPZFNMatrix<9,REAL> kappa_0)
    {
        m_kappa_0 = kappa_0;
    }
    
    /// Get absolute permeability at initial REAL
    TPZFNMatrix<9,REAL> kappa_0()
    {
        return m_kappa_0;
    }
    
    
    /// Set absolute permeability at current REAL
    void Setkappa(TPZFNMatrix<9,REAL> kappa)
    {
        m_kappa = kappa;
    }
    
    /// Get absolute permeability at current REAL
    TPZFNMatrix<9,REAL> kappa()
    {
        return m_kappa;
    }
    
    
    /// Set lagrangian porosity at intial REAL
    void Setphi_0(REAL phi_0)
    {
        m_phi_0 = phi_0;
    }
    
    /// Get lagrangian porosity at intial REAL
    REAL phi_0()
    {
        return m_phi_0;
    }
    
    
    /// Set lagrangian porosity at current REAL
    void Setphi(REAL phi)
    {
        m_phi = phi;
    }
    
    /// Get lagrangian porosity at current REAL
    REAL phi()
    {
        return m_phi;
    }
    
    
    
};


#endif /* TPZMonoPhasicMemoryDFN_h */

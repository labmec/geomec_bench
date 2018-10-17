//
//  TPZMonoPhasicMemoryBCDFN.h
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#ifndef TPZMonoPhasicMemoryBCDFN_h
#define TPZMonoPhasicMemoryBCDFN_h

#include <stdio.h>
#include <iostream>
#include <string>
#include "TPZSavable.h"
#include "Hash/TPZHash.h"
#include "TPZStream.h"
#include "TPZPersistenceManager.h"

class TPZMonoPhasicMemoryBCDFN{
    
private:
    
    /// pore pressure at last state
    STATE m_p;
    
    /// pore pressure at current state
    STATE m_p_n;

    /// normal flux at last state
    STATE m_q;
    
    /// normal flux at current state
    STATE m_q_n;

public:
    
    /// Default constructor
    TPZMonoPhasicMemoryBCDFN();
    
    /// Copy constructor
    TPZMonoPhasicMemoryBCDFN(const TPZMonoPhasicMemoryBCDFN & other);
    
    /// Assignement constructor
    const TPZMonoPhasicMemoryBCDFN & operator=(const TPZMonoPhasicMemoryBCDFN & other);
    
    /// Desconstructor
    virtual ~TPZMonoPhasicMemoryBCDFN();
    
    /// Class name
    const std::string Name() const;
    
    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// Read class attributes
    virtual void Read(TPZStream &buf, void *context);
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    /// Print class attributes
    friend std::ostream & operator<<( std::ostream& out, const TPZMonoPhasicMemoryBCDFN & memory ){
        memory.Print(out);
        return out;
    }
    
    virtual int ClassId() const;
    
    
    /// set and get methods
    
    /// Set pore pressure at current state
    void Setp(STATE p)
    {
        m_p = p;
    }
    
    /// Get pore pressure at current state
    STATE p()
    {
        return m_p;
    }
    
    
    /// Set pore pressure at last state
    void Setp_n(STATE p_n)
    {
        m_p_n = p_n;
    }
    
    /// Get pore pressure at last state
    STATE p_n()
    {
        return m_p_n;
    }

    /// Set normal flux at current state
    void Setq(STATE q)
    {
        m_q = q;
    }
    
    /// Get normal flux at current state
    STATE q()
    {
        return m_q;
    }
    
    
    /// Set normal flux at last state
    void Setq_n(STATE q_n)
    {
        m_q_n = q_n;
    }
    
    /// Get normal flux at last state
    STATE q_n()
    {
        return m_q_n;
    }
    
    
};


#endif /* TPZMonoPhasicMemoryBCDFN_h */

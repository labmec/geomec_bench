//
//  TPZMonoPhasicMemoryDFN.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#include "TPZMonoPhasicMemoryDFN.h"


TPZMonoPhasicMemoryDFN::TPZMonoPhasicMemoryDFN(){
    
    m_p         = 0.0;
    m_p_n       = 0.0;
    m_kappa_0.Resize(3,3);
    m_kappa_0.Zero();
    m_kappa.Resize(3,3);
    m_kappa.Zero();
    m_phi_0     = 0.;
    m_phi       = 0.;
    
}

TPZMonoPhasicMemoryDFN::TPZMonoPhasicMemoryDFN(const TPZMonoPhasicMemoryDFN & other){
    
    m_p         = other.m_p;
    m_p_n       = other.m_p_n;
    m_kappa_0   = other.m_kappa_0;
    m_kappa     = other.m_kappa;
    m_phi_0     = other.m_phi_0;
    m_phi       = other.m_phi;
}

const TPZMonoPhasicMemoryDFN & TPZMonoPhasicMemoryDFN::operator=(const TPZMonoPhasicMemoryDFN & other){
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_p         = other.m_p;
    m_p_n       = other.m_p_n;
    m_kappa_0   = other.m_kappa_0;
    m_kappa     = other.m_kappa;
    m_phi_0     = other.m_phi_0;
    m_phi       = other.m_phi;
    
    return *this;
}

TPZMonoPhasicMemoryDFN::~TPZMonoPhasicMemoryDFN(){
    
}

const std::string TPZMonoPhasicMemoryDFN::Name() const{
    return "TPZMonoPhasicMemoryDFN";
}

void TPZMonoPhasicMemoryDFN::Write(TPZStream &buf, int withclassid) const {
    
    buf.Write(&m_p);
    buf.Write(&m_p_n);
    buf.Write(&m_kappa_0);
    buf.Write(&m_kappa);
    buf.Write(&m_phi_0);
    buf.Write(&m_phi);
}


void TPZMonoPhasicMemoryDFN::Read(TPZStream &buf, void *context){
    
    buf.Read(&m_p);
    buf.Read(&m_p_n);
    buf.Read(&m_phi_0);
    buf.Read(&m_phi);

    for (int i = 0; i< m_kappa_0.Rows(); i++) {
        buf.Read(&m_kappa_0(i,i));
    }
    for (int i = 0; i < m_kappa.Rows(); i++) {
        buf.Read(&m_kappa(i,i));
    }

}

void TPZMonoPhasicMemoryDFN::Print(std::ostream &out) const{
    out << Name();
    out << "\n pressure at last state = " << m_p;
    out << "\n pressure at current state = " << m_p_n;
    out << "\n initial absolute permeability = " << m_kappa_0;
    out << "\n current absolute permeability = " << m_kappa;
    out << "\n initial porosity = " << m_phi_0;
    out << "\n current porosity = " << m_phi;
    out << "\n ";
}

int TPZMonoPhasicMemoryDFN::ClassId() const{
    return Hash("TPZMonoPhasicMemoryDFN");
}

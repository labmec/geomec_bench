//
//  TPZMonoPhasicMemoryBCDFN.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#include "TPZMonoPhasicMemoryBCDFN.h"


TPZMonoPhasicMemoryBCDFN::TPZMonoPhasicMemoryBCDFN(){
    
    m_p         = 0.0;
    m_p_n       = 0.0;
    m_q         = 0.0;
    m_q_n       = 0.0;
}

TPZMonoPhasicMemoryBCDFN::TPZMonoPhasicMemoryBCDFN(const TPZMonoPhasicMemoryBCDFN & other){
    
    m_p         = other.m_p;
    m_p_n       = other.m_p_n;
    m_q         = other.m_q;
    m_q_n       = other.m_q_n;
}

const TPZMonoPhasicMemoryBCDFN & TPZMonoPhasicMemoryBCDFN::operator=(const TPZMonoPhasicMemoryBCDFN & other){
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_p         = other.m_p;
    m_p_n       = other.m_p_n;
    m_q         = other.m_q;
    m_q_n       = other.m_q_n;
    
    return *this;
}

TPZMonoPhasicMemoryBCDFN::~TPZMonoPhasicMemoryBCDFN(){
    
}

const std::string TPZMonoPhasicMemoryBCDFN::Name() const{
    return "TPZMonoPhasicMemoryBCDFN";
}

void TPZMonoPhasicMemoryBCDFN::Write(TPZStream &buf, int withclassid) const {
    buf.Write(&m_p);
    buf.Write(&m_p_n);
    buf.Write(&m_q);
    buf.Write(&m_q_n);
}


void TPZMonoPhasicMemoryBCDFN::Read(TPZStream &buf, void *context){
    buf.Read(&m_p);
    buf.Read(&m_p_n);
    buf.Read(&m_q);
    buf.Read(&m_q_n);
}

void TPZMonoPhasicMemoryBCDFN::Print(std::ostream &out) const{
    out << Name();
    out << "\n pressure at last state = " << m_p;
    out << "\n pressure at current state = " << m_p_n;
    out << "\n normal flux at last state = " << m_p;
    out << "\n normal flux at current state = " << m_p_n;
    out << "\n ";
}

int TPZMonoPhasicMemoryBCDFN::ClassId() const{
    return Hash("TPZMonoPhasicMemoryBCDFN");
}

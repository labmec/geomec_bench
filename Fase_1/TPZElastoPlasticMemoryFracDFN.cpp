//
//  TPZElastoPlasticMemoryFracDFN.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#include "TPZElastoPlasticMemoryFracDFN.h"


TPZElastoPlasticMemoryFracDFN::TPZElastoPlasticMemoryFracDFN(){
    
    m_forceFrac_n.resize(3);
    m_forceFrac_0.resize(3);
    m_u_n.resize(3);
    m_u.resize(3);

    for (int i = 0; i< 3; i++) {
        m_forceFrac_n[i] = 0.;
        m_forceFrac_0[i] = 0.;
        m_u_n[i] = 0.;
        m_u[i] = 0.;
    }
    
}

TPZElastoPlasticMemoryFracDFN::TPZElastoPlasticMemoryFracDFN(const TPZElastoPlasticMemoryFracDFN & other){
   
    m_forceFrac_n       = other.m_forceFrac_n;
    m_forceFrac_0       = other.m_forceFrac_0;
    m_u_n               = other.m_u_n;
    m_u                 = other.m_u;
    
}

const TPZElastoPlasticMemoryFracDFN & TPZElastoPlasticMemoryFracDFN::operator=(const TPZElastoPlasticMemoryFracDFN & other){
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_forceFrac_n       = other.m_forceFrac_n;
    m_forceFrac_0       = other.m_forceFrac_0;
    m_u_n               = other.m_u_n;
    m_u                 = other.m_u;

    return *this;
}

TPZElastoPlasticMemoryFracDFN::~TPZElastoPlasticMemoryFracDFN(){
    
}

const std::string TPZElastoPlasticMemoryFracDFN::Name() const{
    return "TPZElastoPlasticMemoryFracDFN";
}

void TPZElastoPlasticMemoryFracDFN::Write(TPZStream &buf, int withclassid) const {
    buf.Write(m_forceFrac_n);
    buf.Write(m_forceFrac_0);
    buf.Write(m_u_n);
    buf.Write(m_u);
}


void TPZElastoPlasticMemoryFracDFN::Read(TPZStream &buf, void *context){
    buf.Read(m_forceFrac_n);
    buf.Read(m_forceFrac_0);
    buf.Read(m_u_n);
    buf.Read(m_u);
}

void TPZElastoPlasticMemoryFracDFN::Print(std::ostream &out) const{
    out << Name();
    out << "\n Current normal stress in fracture field = " << m_forceFrac_n;
    out << "\n Last normal stress in fracture field = " << m_forceFrac_0;
    out << "\n Current displacement field = " << m_u_n;
    out << "\n Last displacement field = " << m_u;
}

int TPZElastoPlasticMemoryFracDFN::ClassId() const{
    return Hash("TPZElastoPlasticMemoryFracDFN");
}

//
//  TPZElastoPlasticMemoryBCDFN.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#include "TPZElastoPlasticMemoryBCDFN.h"


TPZElastoPlasticMemoryBCDFN::TPZElastoPlasticMemoryBCDFN(){

    m_u.resize(3);
    m_u_n.resize(3);
    m_sigma.Zero();
    m_sigma_n.Zero();
    
}

TPZElastoPlasticMemoryBCDFN::TPZElastoPlasticMemoryBCDFN(const TPZElastoPlasticMemoryBCDFN & other){

    m_u                 = other.m_u;
    m_u_n               = other.m_u_n;
    m_sigma             = other.m_sigma;
    m_sigma_n           = other.m_sigma_n;

}

const TPZElastoPlasticMemoryBCDFN & TPZElastoPlasticMemoryBCDFN::operator=(const TPZElastoPlasticMemoryBCDFN & other){
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }

    m_u                 = other.m_u;
    m_u_n               = other.m_u_n;
    m_sigma             = other.m_sigma;
    m_sigma_n           = other.m_sigma_n;

    return *this;
}

TPZElastoPlasticMemoryBCDFN::~TPZElastoPlasticMemoryBCDFN(){
    
}

const std::string TPZElastoPlasticMemoryBCDFN::Name() const{
    return "TPZElastoPlasticMemoryBCDFN";
}

void TPZElastoPlasticMemoryBCDFN::Write(TPZStream &buf, int withclassid) const {
    buf.Write(m_u);
    buf.Write(m_u_n);
    m_sigma.Write(buf, withclassid);
    m_sigma_n.Write(buf, withclassid);
}


void TPZElastoPlasticMemoryBCDFN::Read(TPZStream &buf, void *context){
    buf.Read(m_u);
    buf.Read(m_u_n);
    m_sigma.Read(buf, context);
    m_sigma_n.Read(buf, context);
}

void TPZElastoPlasticMemoryBCDFN::Print(std::ostream &out) const{
    out << Name();
    out << "\n Last displacement field = " << m_u;
    out << "\n Current displacement field = " << m_u_n;
    out << "\n Last state stress = \n" << m_sigma;
    out << "\n Current state stress = \n" << m_sigma_n;
}

int TPZElastoPlasticMemoryBCDFN::ClassId() const{
    return Hash("TPZElastoPlasticMemoryBCDFN");
}

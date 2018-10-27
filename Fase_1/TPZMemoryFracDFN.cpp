//
//  TPZMemoryFracDFN.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#include "TPZMemoryFracDFN.h"


TPZMemoryFracDFN::TPZMemoryFracDFN() : TPZMonoPhasicMemoryFracDFN() , TPZElastoPlasticMemoryFracDFN() {
    m_alpha = 0.5;
    m_uR.resize(3);
    for (int i=0; i<3; i++) {
        m_uR[i]=0.;
    }
}

TPZMemoryFracDFN::TPZMemoryFracDFN(const TPZMemoryFracDFN & other): TPZMonoPhasicMemoryFracDFN(other), TPZElastoPlasticMemoryFracDFN(other) {
    m_alpha = other.m_alpha;
    m_uR   = other.m_uR;
}

const TPZMemoryFracDFN & TPZMemoryFracDFN::operator=(const TPZMemoryFracDFN & other) {
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_alpha = other.m_alpha;
    m_uR   = other.m_uR;
    
    return *this;
}

TPZMemoryFracDFN::~TPZMemoryFracDFN(){
    
}

const std::string TPZMemoryFracDFN::Name() const {
    return "TPZMemoryFracDFN";
}

void TPZMemoryFracDFN::Write(TPZStream &buf, int withclassid) const {
    TPZMonoPhasicMemoryFracDFN::Write(buf, withclassid);
    TPZElastoPlasticMemoryFracDFN::Write(buf, withclassid);

}

void TPZMemoryFracDFN::Read(TPZStream &buf, void *context){
    TPZMonoPhasicMemoryFracDFN::Read(buf, context);
    TPZElastoPlasticMemoryFracDFN::Read(buf, context);
}

void TPZMemoryFracDFN::Print(std::ostream &out) const {
    out << Name();
    out << "\n Coord of integratrion point = " << m_coord;
    out << "\n";
    TPZMonoPhasicMemoryFracDFN::Print(out);
    TPZElastoPlasticMemoryFracDFN::Print(out);
}

int TPZMemoryFracDFN::ClassId() const {
    return Hash("TPZMemoryFracDFN");
}

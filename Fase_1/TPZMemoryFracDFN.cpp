//
//  TPZMemoryFracDFN.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#include "TPZMemoryFracDFN.h"


TPZMemoryFracDFN::TPZMemoryFracDFN() : TPZMonoPhasicMemoryDFN() , TPZElastoPlasticMemoryDFN() {
    m_alpha = 0.5;
    m_uR.resize(3);

}

TPZMemoryFracDFN::TPZMemoryFracDFN(const TPZMemoryFracDFN & other): TPZMonoPhasicMemoryDFN(other), TPZElastoPlasticMemoryDFN(other) {
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
    TPZMonoPhasicMemoryDFN::Write(buf, withclassid);
    TPZElastoPlasticMemoryDFN::Write(buf, withclassid);

}

void TPZMemoryFracDFN::Read(TPZStream &buf, void *context){
    TPZMonoPhasicMemoryDFN::Read(buf, context);
    TPZElastoPlasticMemoryDFN::Read(buf, context);
}

void TPZMemoryFracDFN::Print(std::ostream &out) const {
    TPZMonoPhasicMemoryDFN::Print(out);
    TPZElastoPlasticMemoryDFN::Print(out);
}

int TPZMemoryFracDFN::ClassId() const {
    return Hash("TPZMemoryFracDFN");
}

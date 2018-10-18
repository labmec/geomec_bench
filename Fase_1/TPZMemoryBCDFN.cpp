//
//  TPZMemoryBCDFN.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#include "TPZMemoryBCDFN.h"


TPZMemoryBCDFN::TPZMemoryBCDFN() : TPZMonoPhasicMemoryBCDFN() , TPZElastoPlasticMemoryBCDFN() {
    m_alpha = 0.5;
}

TPZMemoryBCDFN::TPZMemoryBCDFN(const TPZMemoryBCDFN & other): TPZMonoPhasicMemoryBCDFN(other), TPZElastoPlasticMemoryBCDFN(other) {
    m_alpha = other.m_alpha;
}

const TPZMemoryBCDFN & TPZMemoryBCDFN::operator=(const TPZMemoryBCDFN & other) {
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_alpha = other.m_alpha;

    
    return *this;
}

TPZMemoryBCDFN::~TPZMemoryBCDFN(){
    
}

const std::string TPZMemoryBCDFN::Name() const {
    return "TPZMemoryBCDFN";
}

void TPZMemoryBCDFN::Write(TPZStream &buf, int withclassid) const {
    TPZMonoPhasicMemoryBCDFN::Write(buf, withclassid);
    TPZElastoPlasticMemoryBCDFN::Write(buf, withclassid);
}

void TPZMemoryBCDFN::Read(TPZStream &buf, void *context){
    TPZMonoPhasicMemoryBCDFN::Read(buf, context);
    TPZElastoPlasticMemoryBCDFN::Read(buf, context);
}

void TPZMemoryBCDFN::Print(std::ostream &out) const {
    TPZMonoPhasicMemoryBCDFN::Print(out);
    TPZElastoPlasticMemoryBCDFN::Print(out);
}

int TPZMemoryBCDFN::ClassId() const {
    return Hash("TPZMemoryBCDFN");
}

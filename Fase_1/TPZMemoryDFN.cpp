//
//  TPZMemoryDFN.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#include "TPZMemoryDFN.h"


TPZMemoryDFN::TPZMemoryDFN() : TPZMonoPhasicMemoryDFN() , TPZElastoPlasticMemoryDFN() {
    m_alpha = 0.5;
    m_coord.resize(3);
    
    for (int i = 0; i < m_coord.size(); i++) {
        m_coord[i] = 0.;
    }
}

TPZMemoryDFN::TPZMemoryDFN(const TPZMemoryDFN & other): TPZMonoPhasicMemoryDFN(other), TPZElastoPlasticMemoryDFN(other) {
    m_alpha = other.m_alpha;
    m_coord = other.m_coord;
}

const TPZMemoryDFN & TPZMemoryDFN::operator=(const TPZMemoryDFN & other) {
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_alpha = other.m_alpha;
    m_coord = other.m_coord;

    return *this;
}

TPZMemoryDFN::~TPZMemoryDFN(){
    
}

const std::string TPZMemoryDFN::Name() const {
    return "TPZMemoryDFN";
}

void TPZMemoryDFN::Write(TPZStream &buf, int withclassid) const {
    TPZMonoPhasicMemoryDFN::Write(buf, withclassid);
    TPZElastoPlasticMemoryDFN::Write(buf, withclassid);
}

void TPZMemoryDFN::Read(TPZStream &buf, void *context){
    TPZMonoPhasicMemoryDFN::Read(buf, context);
    TPZElastoPlasticMemoryDFN::Read(buf, context);
}

void TPZMemoryDFN::Print(std::ostream &out) const {
    out << "\n";
    out << "\n -------------------------------";
    out << Name();
    out << "\n Coord of integratrion point = " << m_coord;
    out << "\n -------------------------------";
    
    TPZMonoPhasicMemoryDFN::Print(out);
    out << "\n -------------------------------";
    TPZElastoPlasticMemoryDFN::Print(out);
    out << "\n -------------------------------";

}

int TPZMemoryDFN::ClassId() const {
    return Hash("TPZMemoryDFN");
}

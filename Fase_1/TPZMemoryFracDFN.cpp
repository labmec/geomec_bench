//
//  TPZMemoryFracDFN.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#include "TPZMemoryFracDFN.h"


TPZMemoryFracDFN::TPZMemoryFracDFN() : TPZMonoPhasicMemoryFracDFN() , TPZElastoPlasticMemoryFracDFN() {
    m_alpha = 0.5;
    m_Du_0 = 0.;
    m_Du = 0.;
    m_Du_n = 0.;
    m_Vm = 0.;
    m_coord.resize(3);
    
    for (int i = 0; i < m_coord.size(); i++) {
        m_coord[i] = 0.;
    }

}

TPZMemoryFracDFN::TPZMemoryFracDFN(const TPZMemoryFracDFN & other): TPZMonoPhasicMemoryFracDFN(other), TPZElastoPlasticMemoryFracDFN(other) {
    m_alpha = other.m_alpha;
    m_Du_0 = other.m_Du_0;
    m_Du = other.m_Du;
    m_Du_n = other.m_Du_n;
    m_Vm = other.m_Vm;
    m_coord = other.m_coord;
}

const TPZMemoryFracDFN & TPZMemoryFracDFN::operator=(const TPZMemoryFracDFN & other) {
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_alpha = other.m_alpha;
    m_Du_0 = other.m_Du_0;
    m_Du = other.m_Du;
    m_Du_n = other.m_Du_n;
    m_Vm = other.m_Vm;
    m_coord = other.m_coord;
    
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
    out << "\n";
    out << "\n -------------------------------";
    out << Name();
    out << "\n Coord of integratrion point = " << m_coord;
    out << "\n -------------------------------";
    out << "\n Initial fracture closure = " << m_Du_0;
    out << "\n Last fracture closure = " << m_Du;
    out << "\n Current fracture closure = " << m_Du_n;
    out << "\n Fracture overture = " << m_Vm + m_Du_0 - m_Du_n;
    out << "\n -------------------------------";
    TPZMonoPhasicMemoryFracDFN::Print(out);
    TPZElastoPlasticMemoryFracDFN::Print(out);
}

int TPZMemoryFracDFN::ClassId() const {
    return Hash("TPZMemoryFracDFN");
}

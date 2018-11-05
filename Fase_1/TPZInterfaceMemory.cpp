#include "pzlog.h"
#include "TPZInterfaceMemory.h"

#include <iostream>
#include <cmath>


TPZInterfaceMemory::TPZInterfaceMemory()
{

  m_forceFrac_n.Resize(3);
  m_u_n.Resize(3);
    for (int i=0; i<3; i++) {
        m_forceFrac_n[i] = 0.;
        m_u_n[i] = 0.;
    }

  //this->SetCurrentState();
}


TPZInterfaceMemory::~TPZInterfaceMemory()
{
  
}

/** @brief copy constructor */
TPZInterfaceMemory::TPZInterfaceMemory(const TPZInterfaceMemory &copy) : m_forceFrac_n(copy.m_forceFrac_n), m_u_n(copy.m_u_n)
{
    
}

/** @brief operator equal */
TPZInterfaceMemory & TPZInterfaceMemory::operator=(const TPZInterfaceMemory &copy) {
    m_forceFrac_n = copy.m_forceFrac_n;
    m_u_n = copy.m_u_n;
    return *this;
}

const std::string TPZInterfaceMemory::Name() const{
    return "TPZInterfaceMemory";
}

void TPZInterfaceMemory::Write(TPZStream &buf, int withclassid) const {
    buf.Write(m_forceFrac_n);
    buf.Write(m_u_n);

}


void TPZInterfaceMemory::Read(TPZStream &buf, void *context){
    buf.Read(m_forceFrac_n);
    buf.Read(m_u_n);
}

void TPZInterfaceMemory::Print(std::ostream &out) const{
    out << Name();
    out << "\n current solution at left, forceFrac_n = " << m_forceFrac_n;
    out << "\n current solution at right, u_n = " << m_u_n;
    out << "\n ";
}

int TPZInterfaceMemory::ClassId() const{
    return Hash("TPZInterfaceMemory");
}


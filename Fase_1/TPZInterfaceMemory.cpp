#include "pzlog.h"
#include "TPZInterfaceMemory.h"

#include <iostream>
#include <cmath>


TPZInterfaceMemory::TPZInterfaceMemory() : TPZMonoPhasicMemoryDFN()
{

  SolL.Resize(3);
  SolR.Resize(3);
  //this->SetCurrentState();
}


TPZInterfaceMemory::~TPZInterfaceMemory()
{
  
}

/** @brief copy constructor */
TPZInterfaceMemory::TPZInterfaceMemory(const TPZInterfaceMemory &copy) : TPZMonoPhasicMemoryDFN(copy) {
    
}

/** @brief operator equal */
TPZInterfaceMemory & TPZInterfaceMemory::operator=(const TPZInterfaceMemory &copy) {
    SolL = copy.SolL;
    SolR = copy.SolR;
    return *this;
}

const std::string TPZInterfaceMemory::Name() const{
    return "TPZInterfaceMemory";
}

void TPZInterfaceMemory::Write(TPZStream &buf, int withclassid) const {
    buf.Write(SolL);
    buf.Write(SolR);

}


void TPZInterfaceMemory::Read(TPZStream &buf, void *context){
    buf.Read(SolL);
    buf.Read(SolR);
}

void TPZInterfaceMemory::Print(std::ostream &out) const{
    out << Name();
    out << "\n current solution at left = " << SolL;
    out << "\n current solution at right = " << SolR;
    out << "\n ";
}

int TPZInterfaceMemory::ClassId() const{
    return Hash("TPZInterfaceMemory");
}


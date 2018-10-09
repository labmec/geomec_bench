//
//  TPZMatFractureBB.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/2/14.
//
//

#include "TPZMatFractureBB.h"
#include "pzaxestools.h"


template <class TMEM>
int TPZMatFractureBB<TMEM>::ClassId() const{
    return Hash("TPZMatFractureBB") ^ TPZDiscontinuousGalerkin::ClassId() << 1;
}

template <class TMEM>
void TPZMatFractureBB<TMEM>::Write(TPZStream &buf, int withclassid) const
{
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    buf.Write(&fNStateVariables);
}

template <class TMEM>
void TPZMatFractureBB<TMEM>::Read(TPZStream &buf, void *context)
{
    TPZDiscontinuousGalerkin::Read(buf, context);
    buf.Read(&fNStateVariables);
    
}

template <class TMEM>
void TPZMatFractureBB<TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    ek.Zero();
    ef.Zero();
    
    return;
    DebugStop();
}

template <class TMEM>
void TPZMatFractureBB<TMEM>::UpdateMemory(TPZVec<TPZMaterialData> &datavec){
    const int intGlobPtIndex = datavec[0].intLocPtIndex;
//    TPZFMatrix<REAL> Vl = this->MemItem(intGlobPtIndex);
//    const STATE pfrac = datavec[1].sol[0][0];
//    const REAL deltaT = fData->TimeStep();
//    REAL tStar = fData->FictitiousTime(Vl(0,0), pfrac);
//    REAL Vlnext = fData->VlFtau(pfrac, tStar + deltaT);
//    Vl(0,0) = Vlnext;
//    this->MemItem(intGlobPtIndex) = Vl;
}

template <class TMEM>
void TPZMatFractureBB<TMEM>::UpdateMemory(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavec){
    const int intGlobPtIndex = data.intGlobPtIndex;
//    TPZFMatrix<REAL> Vl = this->MemItem(intGlobPtIndex);
//    const STATE pfrac = datavec[1].sol[0][0];
//    const REAL deltaT = fData->TimeStep();
//    REAL tStar = fData->FictitiousTime(Vl(0,0), pfrac);
//    REAL Vlnext = fData->VlFtau(pfrac, tStar + deltaT);
//    Vl(0,0) = Vlnext;
//    this->MemItem(intGlobPtIndex) = Vl;
}



template class TPZMatFractureBB<TPZMemoryFracDFN>;

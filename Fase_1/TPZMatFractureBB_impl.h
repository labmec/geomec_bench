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
void TPZMatFractureBB<TMEM>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    int nmesh = datavec.size();
    if (nmesh!=2) DebugStop();
    
    TPZFMatrix<REAL>  &phiQ = datavec[0].phi;
    TPZFMatrix<REAL> &phiP = datavec[1].phi;
    int phrq = phiQ.Rows();
    int phrp = phiP.Rows();
    
    //------- Block of matrix B ------
    int iq, jp;
    for(iq = 0; iq<phrq; iq++) {
        for(jp=0; jp<phrp; jp++) {
            ek(iq, phrq+jp) += fMultiplier*weight*phiQ(iq,0)*phiP(jp,0);
        }
    }
    
    
    //------- Block of matrix B^T ------
    int ip, jq;
    for(ip=0; ip<phrp; ip++) {
        for(jq=0; jq<phrq; jq++) {
            ek(ip + phrq,jq) += fMultiplier*weight*phiP(ip,0)*phiQ(jq,0);
        }
    }
}

template <class TMEM>
void TPZMatFractureBB<TMEM>::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    TPZFMatrix<REAL> *phiLPtr = 0, *phiRPtr = 0;
    for (int i=0; i<dataleft.size(); i++) {
        if (dataleft[i].phi.Rows() != 0) {
            phiLPtr = &dataleft[i].phi;
            break;
        }
    }
    for (int i=0; i<dataright.size(); i++) {
        if (dataright[i].phi.Rows() != 0) {
            phiRPtr = &dataright[i].phi;
            break;
        }
    }
    
    if(!phiLPtr || !phiRPtr)
    {
        DebugStop();
    }
    TPZFMatrix<REAL> &phiL = *phiLPtr;
    TPZFMatrix<REAL> &phiR = *phiRPtr;
    
    
    int nrowl = phiL.Rows();
    int nrowr = phiR.Rows();
    
    if(nrowl+nrowr != ek.Rows())
    {
        std::cout<<ek.Rows()<<std::endl;
        std::cout<<nrowl<<std::endl;
        std::cout<<nrowr<<std::endl;
        DebugStop();
    }
    
    int secondblock = ek.Rows()-phiR.Rows()*fNStateVariables;
    int il,jl,ir,jr;
    
    // 3) phi_I_left, phi_J_right
    for(il=0; il<nrowl; il++) {
        for(jr=0; jr<nrowr; jr++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ek(fNStateVariables*il+ist,fNStateVariables*jr+ist+secondblock) += weight * fMultiplier * (phiL(il) * phiR(jr));
            }
        }
    }
    
    //    // 4) phi_I_right, phi_J_left
    for(ir=0; ir<nrowr; ir++) {
        for(jl=0; jl<nrowl; jl++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ek(ir*fNStateVariables+ist+secondblock,jl*fNStateVariables+ist) += weight * fMultiplier * (phiR(ir) * phiL(jl));
            }
        }
    }
    
}


template <class TMEM>
void TPZMatFractureBB<TMEM>::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    //    TPZFMatrix<REAL> &dphiLdAxes = dataleft.dphix;
    //    TPZFMatrix<REAL> &dphiRdAxes = dataright.dphix;
    TPZFMatrix<REAL> &phiL = dataleft.phi;
    TPZFMatrix<REAL> &phiR = dataright.phi;
    
    //    TPZFNMatrix<660> dphiL, dphiR;
    //    TPZAxesTools<REAL>::Axes2XYZ(dphiLdAxes, dphiL, dataleft.axes);
    //    TPZAxesTools<REAL>::Axes2XYZ(dphiRdAxes, dphiR, dataright.axes);
    
    int nrowl = phiL.Rows();
    int nrowr = phiR.Rows();
#ifdef PZDEBUG
    if(phiL.Rows()*fNStateVariables+phiR.Rows()*fNStateVariables != ek.Rows())
    {
        DebugStop();
    }
#endif
    int secondblock = ek.Rows()-phiR.Rows()*fNStateVariables;
    int il,jl,ir,jr;
    
    // 3) phi_I_left, phi_J_right
    for(il=0; il<nrowl; il++) {
        for(jr=0; jr<nrowr; jr++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ek(fNStateVariables*il+ist,fNStateVariables*jr+ist+secondblock) += weight * fMultiplier * (phiL(il) * phiR(jr));
            }
        }
    }
    
    //    // 4) phi_I_right, phi_J_left
    for(ir=0; ir<nrowr; ir++) {
        for(jl=0; jl<nrowl; jl++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ek(ir*fNStateVariables+ist+secondblock,jl*fNStateVariables+ist) += weight * fMultiplier * (phiR(ir) * phiL(jl));
            }
        }
    }
    
}

template <class TMEM>
void TPZMatFractureBB<TMEM>::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef)
{
    
    TPZFMatrix<REAL> &phiL = dataleft.phi;
    TPZFMatrix<REAL> &phiR = dataright.phi;
    
    int nrowl = phiL.Rows();
    int nrowr = phiR.Rows();
    
    STATE vn    = dataleft.sol[0][0];
    int secondblock = ef.Rows()-phiR.Rows()*fNStateVariables;
 
    // 3) phi_I_left, phi_J_right
    for(int il=0; il<nrowl; il++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ef(fNStateVariables*il+ist) += weight * fMultiplier * (phiL(il) * vn);
            }
    }


    for(int ir=0; ir<nrowr; ir++) {
        for (int ist=0; ist<fNStateVariables; ist++) {
            ef(fNStateVariables*ir+ist+secondblock) += weight * fMultiplier * (phiR(ir) * vn);
        }
    }
    
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

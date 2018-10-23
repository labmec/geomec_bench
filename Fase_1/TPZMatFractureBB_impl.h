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
void TPZMatFractureBB<TMEM>::Write(TPZStream &buf, int withclassid)
{
    TPZMaterial::Write(buf, withclassid);
}

template <class TMEM>
void TPZMatFractureBB<TMEM>::Read(TPZStream &buf, void *context)
{
    TPZMaterial::Read(buf, context);
}

template <class TMEM>
void TPZMatFractureBB<TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    //if(data.intGlobPtIndex == -1) DebugStop();
    
   // long gp_index = data.intGlobPtIndex;
   // TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    
//    TPZManVector<STATE,3> vLn2    = memory.GetLeftSol();
//    TPZManVector<STATE,3> vRn2    = memory.GetRightSol();
//
//    TPZManVector<STATE,3> vLn    = dataleft.sol[0];
//    TPZManVector<STATE,3> vRn    = dataright.sol[0];
//
//    for(int i=0; i<fNStateVariables; i++)
//    {
//        vLn[i] += vLn2[i];
//        vRn[i] += vRn2[i];
//    }
//
//    TPZFMatrix<REAL> &phiL = dataleft.phi;
//    TPZFMatrix<REAL> &phiR = dataright.phi;
//
//    //    TPZFNMatrix<660> dphiL, dphiR;
//    //    TPZAxesTools<REAL>::Axes2XYZ(dphiLdAxes, dphiL, dataleft.axes);
//    //    TPZAxesTools<REAL>::Axes2XYZ(dphiRdAxes, dphiR, dataright.axes);
//
//    int nrowl = phiL.Rows();
//    int nrowr = phiR.Rows();
//#ifdef PZDEBUG
//    if(phiL.Rows()*fNStateVariables+phiR.Rows()*fNStateVariables != ek.Rows())
//    {
//        DebugStop();
//    }
//#endif
//    int secondblock = ek.Rows()-phiR.Rows()*fNStateVariables;
//    int il,jl,ir,jr;
//
//
//
//    // 3) phi_I_left, phi_J_right
//    for(il=0; il<nrowl; il++) {
//
//        for (int ist=0; ist<fNStateVariables; ist++) {
//            ef(fNStateVariables*il+ist) += weight * fMultiplier * (phiL(il) * vRn[ist]);
//        }
//
//        for(jr=0; jr<nrowr; jr++) {
//            for (int ist=0; ist<fNStateVariables; ist++) {
//                ek(fNStateVariables*il+ist,fNStateVariables*jr+ist+secondblock) += weight * fMultiplier * (phiL(il) * phiR(jr));
//            }
//        }
//    }
//
//    //    // 4) phi_I_right, phi_J_left
//    for(ir=0; ir<nrowr; ir++) {
//
//        for (int ist=0; ist<fNStateVariables; ist++) {
//            ef(fNStateVariables*ir+ist+secondblock) += weight * fMultiplier * (phiR(ir) * vLn[ist]);
//        }
//
//        for(jl=0; jl<nrowl; jl++) {
//            for (int ist=0; ist<fNStateVariables; ist++) {
//                ek(ir*fNStateVariables+ist+secondblock,jl*fNStateVariables+ist) += weight * fMultiplier * (phiR(ir) * phiL(jl));
//            }
//        }
//    }
//
//
    ek.Zero();
    ef.Zero();
    
    return;
    DebugStop();
}

template <class TMEM>
void TPZMatFractureBB<TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
   
    return;
    
    if (m_simulation_data->Get_must_accept_solution_Q()) {
        long gp_index = data.intGlobPtIndex;
        if(gp_index < 0)
        {
            std::cout << "Integration point index is not initialized\n";
            DebugStop();
        }
        TPZManVector<STATE,3> sol = data.sol[0];

        TMEM &mem = this->GetMemory().get()->operator[](gp_index);
        TPZManVector<STATE,3> forceFrac_n = mem.GetForceFrac_n();
        for(int i=0; i<fNStateVariables; i++)
        {
            forceFrac_n[i] += sol[i];
        }
        mem.SetForceFrac_n(forceFrac_n);

        
    }else{
        TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
        this->Contribute(data, weight, ek_fake, ef);
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

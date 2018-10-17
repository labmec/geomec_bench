//
//  TPZLagrangeInterface.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/2/14.
//
//

#include "TPZLagrangeInterface.h"
#include "pzaxestools.h"


template <class TMEM>
int TPZLagrangeInterface<TMEM>::ClassId() const{
    return Hash("TPZLagrangeInterface") ^ TPZDiscontinuousGalerkin::ClassId() << 1;
}

template <class TMEM>
void TPZLagrangeInterface<TMEM>::Write(TPZStream &buf, int withclassid) const
{
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    buf.Write(&fNStateVariables);
}

template <class TMEM>
void TPZLagrangeInterface<TMEM>::Read(TPZStream &buf, void *context)
{
    TPZDiscontinuousGalerkin::Read(buf, context);
    buf.Read(&fNStateVariables);
    
}

template <class TMEM>
void TPZLagrangeInterface<TMEM>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
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
void TPZLagrangeInterface<TMEM>::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    
    long gp_index = data.intGlobPtIndex;
   // TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    
    
    int leftdataindex = -1;
    TPZFMatrix<REAL> *phiLPtr = 0, *phiRPtr = 0;
    for (int i=0; i<dataleft.size(); i++) {
        if (dataleft[i].phi.Rows() != 0) {
            phiLPtr = &dataleft[i].phi;
            leftdataindex = i;
            break;
        }
    }
    int rightdataindex = -1;
    for (int i=0; i<dataright.size(); i++) {
        if (dataright[i].phi.Rows() != 0) {
            phiRPtr = &dataright[i].phi;
            rightdataindex = i;
            break;
        }
    }
    
    if(!phiLPtr || !phiRPtr)
    {
        DebugStop();
    }
    TPZFMatrix<REAL> &phiL = *phiLPtr;
    TPZFMatrix<REAL> &phiR = *phiRPtr;
    
    TPZManVector<STATE> vLn    = dataleft[leftdataindex].sol[0];
    TPZManVector<STATE> vRn    = dataright[rightdataindex].sol[0];
    
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
        for (int ist=0; ist<fNStateVariables; ist++) {
            ef(fNStateVariables*il+ist) += weight * fMultiplier * (phiL(il) * vRn[ist]);
        }
        for(jr=0; jr<nrowr; jr++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ek(fNStateVariables*il+ist,fNStateVariables*jr+ist+secondblock) += weight * fMultiplier * (phiL(il) * phiR(jr));
            }
        }
    }
    
    //    // 4) phi_I_right, phi_J_left
    for(ir=0; ir<nrowr; ir++) {
        
        for (int ist=0; ist<fNStateVariables; ist++) {
            ef(fNStateVariables*ir+ist+secondblock) += weight * fMultiplier * (phiR(ir) * vLn[ist]);
        }
        
        for(jl=0; jl<nrowl; jl++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ek(ir*fNStateVariables+ist+secondblock,jl*fNStateVariables+ist) += weight * fMultiplier * (phiR(ir) * phiL(jl));
            }
        }
    }

    
}

template <class TMEM>
void TPZLagrangeInterface<TMEM>::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    //    TPZFMatrix<REAL> &dphiLdAxes = dataleft.dphix;
    //    TPZFMatrix<REAL> &dphiRdAxes = dataright.dphix;
    
    long gp_index = data.intGlobPtIndex;
    TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    
    TPZManVector<STATE,3> vLn2    = memory.GetLeftSol();
    TPZManVector<STATE,3> vRn2    = memory.GetRightSol();
    
    TPZManVector<STATE,3> vLn    = dataleft.sol[0];
    TPZManVector<STATE,3> vRn    = dataright.sol[0];

    for(int i=0; i<fNStateVariables; i++)
    {
        vLn[i] += vLn2[i];
        vRn[i] += vRn2[i];
    }
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
        
        for (int ist=0; ist<fNStateVariables; ist++) {
            ef(fNStateVariables*il+ist) += weight * fMultiplier * (phiL(il) * vRn[ist]);
        }
        
        for(jr=0; jr<nrowr; jr++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ek(fNStateVariables*il+ist,fNStateVariables*jr+ist+secondblock) += weight * fMultiplier * (phiL(il) * phiR(jr));
            }
        }
    }
    
    //    // 4) phi_I_right, phi_J_left
    for(ir=0; ir<nrowr; ir++) {
        
        for (int ist=0; ist<fNStateVariables; ist++) {
            ef(fNStateVariables*ir+ist+secondblock) += weight * fMultiplier * (phiR(ir) * vLn[ist]);
        }
        
        for(jl=0; jl<nrowl; jl++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ek(ir*fNStateVariables+ist+secondblock,jl*fNStateVariables+ist) += weight * fMultiplier * (phiR(ir) * phiL(jl));
            }
        }
    }

}

template <class TMEM>
void TPZLagrangeInterface<TMEM>::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef)
{
    
    if (m_simulation_data->Get_must_accept_solution_Q()) {
        long gp_index = data.intGlobPtIndex;
        if(gp_index < 0)
        {
            std::cout << "Integration point index is not initialized\n";
            DebugStop();
        }
        TPZManVector<STATE,3> vL    = dataleft.sol[0];
        TPZManVector<STATE,3> vR    = dataright.sol[0];

//        if (m_simulation_data->IsInitialStateQ()) {
//            this->GetMemory().get()->operator[](gp_index).Setp_0(p);
//            this->GetMemory().get()->operator[](gp_index).Setp_0(p);
//        }

//        if (m_simulation_data->IsCurrentStateQ()) {
//            this->GetMemory().get()->operator[](gp_index).SetRightSol_n(vR);
//            this->GetMemory().get()->operator[](gp_index).SetLeftSol_n(vL);
//        }else{
        TMEM &mem = this->GetMemory().get()->operator[](gp_index);
        TPZManVector<STATE,3> solL = mem.GetLeftSol();
        TPZManVector<STATE,3> solR = mem.GetRightSol();
        for(int i=0; i<fNStateVariables; i++)
        {
            solL[i] += vL[i];
            solR[i] += vR[i];
        }
        mem.SetLeftSol(solL);
        mem.SetRightSol(solR);
//        }

    }else{
        TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
        this->ContributeInterface(data, dataleft, dataright, weight, ek_fake, ef);
    }


    return;
    
    
    TPZFMatrix<REAL> &phiL = dataleft.phi;
    TPZFMatrix<REAL> &phiR = dataright.phi;
    
    int nrowl = phiL.Rows();
    int nrowr = phiR.Rows();
    
    TPZManVector<STATE> vLn    = dataleft.sol[0];
    TPZManVector<STATE> vRn    = dataright.sol[0];
    
    int secondblock = ef.Rows()-phiR.Rows()*fNStateVariables;
 
    // 3) phi_I_left, phi_J_right
    for(int il=0; il<nrowl; il++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ef(fNStateVariables*il+ist) += weight * fMultiplier * (phiL(il) * vRn[ist]);
            }
    }


    for(int ir=0; ir<nrowr; ir++) {
        for (int ist=0; ist<fNStateVariables; ist++) {
            ef(fNStateVariables*ir+ist+secondblock) += weight * fMultiplier * (phiR(ir) * vLn[ist]);
        }
    }
    
}

template <class TMEM>
void TPZLagrangeInterface<TMEM>::UpdateMemory(TPZVec<TPZMaterialData> &datavec){
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
void TPZLagrangeInterface<TMEM>::UpdateMemory(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavec){
    const int intGlobPtIndex = data.intGlobPtIndex;
//    TPZFMatrix<REAL> Vl = this->MemItem(intGlobPtIndex);
//    const STATE pfrac = datavec[1].sol[0][0];
//    const REAL deltaT = fData->TimeStep();
//    REAL tStar = fData->FictitiousTime(Vl(0,0), pfrac);
//    REAL Vlnext = fData->VlFtau(pfrac, tStar + deltaT);
//    Vl(0,0) = Vlnext;
//    this->MemItem(intGlobPtIndex) = Vl;
}



template class TPZLagrangeInterface<TPZInterfaceMemory>;

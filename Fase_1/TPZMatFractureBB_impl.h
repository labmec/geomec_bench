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
    
    STATE D_u0 =0; //Initial closute
    STATE Vm = Get_Vm(); //Max opening
    STATE a0 = Get_a0(); //Initial opening
    STATE Kni = Get_Kni(); //Initial normal stiffness
    D_u0 = Vm-a0;
    
    
    long gp_index = data.intGlobPtIndex;
    
    TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    memory.SetCoord(data.x);
    
    
    std::ofstream fileMemFrac("MemoryFrac.txt");
    memory.Print(fileMemFrac);
    
    
    TPZManVector<STATE,3> delta_forceFrac_n    = memory.GetForceFrac_n();

    TPZManVector<STATE,3> forceFrac_n    = data.sol[0];

    for(int i=0; i<2; i++)
    {
        forceFrac_n[i] += delta_forceFrac_n[i];
    }

    TPZFMatrix<REAL> &phi = data.phi;
    TPZFMatrix<REAL> &dphi = data.dphi;
    
    //    TPZFNMatrix<660> dphiL, dphiR;
    //    TPZAxesTools<REAL>::Axes2XYZ(dphiLdAxes, dphiL, dataleft.axes);
    //    TPZAxesTools<REAL>::Axes2XYZ(dphiRdAxes, dphiR, dataright.axes);

    int nrow = phi.Rows();

#ifdef PZDEBUG
    if(phi.Rows()*fNStateVariables != ek.Rows())
    {
        DebugStop();
    }
#endif
    TPZFNMatrix<3,STATE> phiVi(3,1,0.0),phiVj(3,1,0.0);
    // 3) phi_I_left, phi_J_right
    for(int i=0; i<nrow; i++) {

    //    int iphi = data.fShapeIndex[i].second;
    //    int ivec = data.fVecShapeIndex[i].first;
        

        ef(fNStateVariables*i+1) += weight * (phi(i,0) * 100.);
        for (int e = 0; e < 2; e++) {
            phiVi(e,0) = phi(i,0)*data.normal[e];
//            ef(fNStateVariables*i+e) += -weight * (phiVi(e,0) * D_u0);

            REAL valLinear = (weight * phiVi(e,0) * forceFrac_n[e]);
            REAL val = (weight * phiVi(e,0) * forceFrac_n[e] * Vm)/(forceFrac_n[e]+Kni*Vm);
        
            ef(fNStateVariables*i+e) += valLinear*0.;
            
            TPZFNMatrix<3,STATE> phiVj(3,1,0.0);
            for(int j=0; j<nrow; j++) {
        
                for (int f = 0; f < 2; f++) {
                    phiVj(f,0) = phi(j,0)*data.normal[f];

                    REAL valLinear =  weight * (phiVi(e,0) * phiVj(f,0) );
                    REAL val = weight * (phiVi(e,0) * phiVj(f,0) * Vm)/(phiVj(f,0)+Kni*Vm);
                    
                    ek(fNStateVariables*i+e,fNStateVariables*j+f) += valLinear*0.;
                }
                
            }
        
        }


        
//        ef(fNStateVariables*i+1) += weight * (phi(i,0) * 100.);
        
//        for (int e = 0; e < 2; e++) {
//            phiVi(e,0) = phi(i,0)*data.normal[e];
//            ef(fNStateVariables*i+ist) += -weight * (phiVi(e,0) * D_u0);

//            ef(fNStateVariables*i+ist) += weight * (phiVi(e,0) * forceFrac_n[ist] * Vm)/(forceFrac_n[ist]+Kni*Vm);
//
//            for(int j=0; j<nrow; j++) {
//                phiVj(e,0) = phi(j,0)*data.normal[e];
//            for (int e = 0; e < 2; e++) {
//                ek(fNStateVariables*i+ist,fNStateVariables*j+ist+secondblock) += weight * (phiVi(e,0) * phiVj(e,0) * Vm)/(phi(j)+Kni*Vm);
//            }
//            }
      }


}

template <class TMEM>
void TPZMatFractureBB<TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
   
    //return;
    
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

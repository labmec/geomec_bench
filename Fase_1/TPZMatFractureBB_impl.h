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
    
    STATE D_u0 = 0; //Initial closute
    STATE Vm = Get_Vm(); //Max opening
    STATE a0 = Get_a0(); //Initial opening
    STATE Kni = Get_Kni(); //Initial normal stiffness
    
    long gp_index = data.intGlobPtIndex;
    
    TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    memory.SetCoord(data.x);
    D_u0 = memory.GetDu_0();
    
    
    std::ofstream fileMemFrac("MemoryFrac.txt", std::ofstream::app);
    //memory.Print(fileMemFrac);
    
    TPZManVector<STATE,3> delta_forceFrac_n    = memory.GetForceFrac_n();
    TPZManVector<STATE,3> forceFrac_n    = data.sol[0];
    TPZManVector<STATE,3> normal = data.normal;
    
    for(int i=0; i<2; i++)
    {
        forceFrac_n[i] += delta_forceFrac_n[i];
    }

    fileMemFrac << "Coordenadas = " << data.x << std::endl;
    fileMemFrac << "forceFrac = " << forceFrac_n << std::endl;
    
    STATE forceFrac_normal = InnerVec(forceFrac_n, normal);
    
    TPZFMatrix<REAL> &phi_f = data.phi;
    TPZFMatrix<REAL> &dphi = data.dphi;
    
    //    TPZFNMatrix<660> dphiL, dphiR;
    //    TPZAxesTools<REAL>::Axes2XYZ(dphiLdAxes, dphiL, dataleft.axes);
    //    TPZAxesTools<REAL>::Axes2XYZ(dphiRdAxes, dphiR, dataright.axes);

    int nrow = phi_f.Rows();

#ifdef PZDEBUG
    if(phi_f.Rows()*fNStateVariables != ek.Rows())
    {
        DebugStop();
    }
#endif
    TPZManVector<STATE,3> phiVi(3,0.0),phiVj(3,0.0);
    STATE phiVi_normal = 0., phiVj_normal = 0.;
    
    for(int i=0; i<nrow; i++) {

    // int iphi = data.fShapeIndex[i].second;
    // int ivec = data.fVecShapeIndex[i].first;

        //ef(fNStateVariables*i+1) += weight * (phi(i,0) * 100.);
        for (int e = 0; e < 2; e++) {
            phiVi[e] = phi_f(i,0)*normal[e];
        }
        phiVi_normal = InnerVec(phiVi, normal);
        
        // ef(fNStateVariables*i+1) += -weight * (phiVi_normal * 0.);
        // ef(fNStateVariables*i) += -weight * (phiVi_normal * D_u0);

        for(int ist=0; ist<2; ist++)
        {

            STATE valBB = weight * phi_f(i) * forceFrac_n[ist] * Vm /(forceFrac_n[ist] + Kni*Vm);
            
            STATE valNonLinear = weight * phi_f(i) * (forceFrac_n[ist] * forceFrac_n[ist])/200.;
            STATE valLinear = weight * phi_f(i) * forceFrac_n[ist];
            
            STATE valDu_0 = weight * phi_f(i) * D_u0 * normal[ist];
            
            ef(fNStateVariables*i+ist) += -valBB+valDu_0;
        }
        
        //STATE val = (weight * phiVi_normal * forceFrac_normal * Vm)/(forceFrac_normal+Kni*Vm);
        

            
        TPZManVector<STATE,3> phiVj(3,0.0);
        for(int j=0; j<nrow; j++) {
                
            for (int f = 0; f < 2; f++) {
                phiVj[f] = phi_f(j,0)*normal[f];
            }
            phiVj_normal = InnerVec(phiVj, normal);
            
            for(int ist=0; ist<2; ist++)
            {

                STATE valBB = weight * phi_f(i) * phi_f(j) * Kni * Vm * Vm / ( (Kni*Vm+ forceFrac_n[ist]) * (Kni*Vm+ forceFrac_n[ist]) );
                
                STATE valNonLinear =  weight * phi_f(i) * phi_f(j) * 2. * (forceFrac_n [ist])/200.;
                 //STATE valLinear =  weight * phi_f(i) * phi_f(j);
                 ek(fNStateVariables*i+ist,fNStateVariables*j+ist) += -valBB;
            }
           
            // STATE val = weight * (phiVi_normal * phiVj_normal * Vm)/(phiVj_normal+Kni*Vm);
            //2. * forceFrac_n[ist] * * forceFrac_n[ist]

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
    
  //  ek.Print(cout);
  //  ef.Print(cout);
    
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
        TPZManVector<STATE,3> delta_forceFrac = data.sol[0];
        TPZManVector<STATE,3> normal = data.normal;
        
        STATE forceFrac_normal = InnerVec(delta_forceFrac,normal);
        
        STATE Vm = Get_Vm(); //Max opening
        STATE a0 = Get_a0(); //Initial opening
        STATE Kni = Get_Kni(); //Initial normal stiffness
        
        if (m_simulation_data->IsInitialStateQ()) {
            TMEM &mem = this->GetMemory().get()->operator[](gp_index);
            STATE Du = (forceFrac_normal * Vm)/(forceFrac_normal+Kni*Vm);
            mem.SetDu_0(Du);
        }
        
        STATE Du_0 = Vm - a0;
        STATE Du_n = (forceFrac_normal * Vm)/(forceFrac_normal+Kni*Vm);
        
        TMEM &mem = this->GetMemory().get()->operator[](gp_index);
        
        mem.SetVm(Vm);
        mem.SetDu_n(Du_n);
        
        TPZManVector<STATE,3> forceFrac_n = mem.GetForceFrac_n();
        for(int i=0; i<fNStateVariables; i++)
        {
            forceFrac_n[i] += delta_forceFrac[i];
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


template <class TMEM>
STATE TPZMatFractureBB<TMEM>::InnerVec(TPZManVector<STATE,3>  &S, TPZManVector<STATE,3>  &T){
    
    //inner product of two vectors
    
    STATE Val=0.;
    
    for(int j = 0; j < S.size(); j++){
            Val += S[j]*T[j];
    }
    
    return Val;
    
}


template class TPZMatFractureBB<TPZMemoryFracDFN>;

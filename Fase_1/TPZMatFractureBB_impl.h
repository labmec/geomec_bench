//
//  TPZMatFractureBB.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/2/14.
//
//

#include "TPZMatFractureBB.h"
#include "pzaxestools.h"

#ifdef LOG4CXX
static LoggerPtr BartonBandisLogger(Logger::getLogger("Benchmark.BartonBandis"));
#endif

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

   // STATE Kni = Get_Kni(); //Initial normal stiffness
    
    long gp_index = data.intGlobPtIndex;
    
    TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    memory.SetCoord(data.x);
    STATE  D_u0 = memory.GetDu_0(); //Initial closure;
    //STATE  D_un = memory.GetDu_n();
    REAL Vm = memory.GetVm(); //Max opening
    REAL Kni = m_simulation_data->Get_Kni();
    //Kni = 12041;
    
    std::ofstream fileMemFrac("MemoryFrac.txt", std::ofstream::app);
    //memory.Print(fileMemFrac);
    
   // TPZManVector<STATE,3> delta_forceFrac_n    = memory.GetForceFrac_n();
    TPZManVector<STATE,3> forceFrac_n    = data.sol[0];
    TPZManVector<STATE,3> normal = data.normal;
    Correctnormal(normal);
    
//    for(int i=0; i<2; i++)
//    {
//        forceFrac_n[i] += delta_forceFrac_n[i];
//    }

    fileMemFrac << "Coordenadas = " << data.x << std::endl;
    fileMemFrac << "forceFrac = " << forceFrac_n << std::endl;

    
    STATE forceFrac_normal = InnerVec(forceFrac_n, normal);
    STATE p_n = memory.p_n();
    STATE forceFrac_normal_Ef = forceFrac_normal - p_n; //alpha multiplica p
    
    for (int i =0 ; i<forceFrac_n.size(); i++) {
        forceFrac_n[i] = forceFrac_normal_Ef * normal[i];
    }
    
    
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
        
        for (int e = 0; e < 2; e++) {
            phiVi[e] = phi_f(i,0)*normal[e];
        }
        phiVi_normal = InnerVec(phiVi, normal);
        

        for(int ist=0; ist<2; ist++)
        {

            STATE valBB = weight * phi_f(i,0)*normal[ist] * forceFrac_normal_Ef * Vm /(forceFrac_normal_Ef + Kni*Vm);
            //valBB =0.;
            STATE valDu_0 = weight * phi_f(i,0)*normal[ist]  * D_u0 ;
            //valDu_0 = 0.;
            
            ef(fNStateVariables*i+ist) += -valBB+valDu_0;
        }
        

            
        TPZManVector<STATE,3> phiVj(3,0.0);
        for(int j=0; j<nrow; j++) {
                
            for (int f = 0; f < 2; f++) {
                phiVj[f] = phi_f(j,0)*normal[f];
            }
            phiVj_normal = InnerVec(phiVj, normal);
            
                for(int ist=0; ist<2; ist++)
                {

                    for(int jst=0;jst<2;jst++){
                        STATE valBB = weight * phi_f(i,0)*normal[ist] * phi_f(j,0)*normal[jst] * Kni * Vm * Vm / ( (Kni*Vm+ forceFrac_normal_Ef) * (Kni*Vm+ forceFrac_normal_Ef) );
                        //valBB =0.;
                        ek(fNStateVariables*i+ist,fNStateVariables*j+ist) += -valBB;
                    }

                    
//                    if(valBB!=valBB){
//                        DebugStop();
//                    }
                }
           

            }
        
        }

    
#ifdef LOG4CXX
    if (BartonBandisLogger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "<<< TPZMatFractureBB<TMEM>::Contribute ***";
        sout << " Resultant rhs vector:\n" << ef;
        sout << " Resultant stiff vector:\n" << ek;
        LOGPZ_DEBUG(BartonBandisLogger, sout.str().c_str());
    }
#endif

    
}

template <class TMEM>
void TPZMatFractureBB<TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
   
    
    
    if (m_simulation_data->Get_must_accept_solution_Q()) {
        long gp_index = data.intGlobPtIndex;
        
        
        if(gp_index < 0)
        {
            std::cout << "Integration point index is not initialized\n";
            
            DebugStop();
        }
        TPZManVector<STATE,3> forceFrac = data.sol[0];
        TPZManVector<STATE,3> normal = data.normal;
        Correctnormal(normal);
        
        STATE forceFrac_normal = InnerVec(forceFrac,normal);

        TMEM &mem = this->GetMemory().get()->operator[](gp_index);
        STATE  D_u0 = mem.GetDu_0(); //Initial closure;
        REAL Vm = mem.GetVm(); //Max opening
        REAL Kni = m_simulation_data->Get_Kni();

        STATE p_n = mem.p_n();
        STATE forceFrac_normal_Ef = forceFrac_normal - p_n;
        int frac_Id = Frac_ID();
        if (m_simulation_data->IsInitialStateQ()) {
            
            
           // mem.sigma_0();
            int size = m_simulation_data->Get_Stress0().Rows();
            TPZManVector<STATE,3> forceFrac0(3,0.);
            for(int i = 0 ; i<size; i++){
                forceFrac0[i]=m_simulation_data->Get_Stress0()(i,i)*normal[i];
            }
            STATE forceFrac_normal0 = InnerVec(forceFrac0,normal);
            //forceFrac_normal0 = -128.41549701748323;
            STATE Du_0 = (forceFrac_normal0 * Vm)/(forceFrac_normal0+Kni*Vm); //Valor da tensão imposta
           // Du_0 = 0.;

            mem.SetForceFrac_normal_0(forceFrac_normal0);
            mem.SetDu_0(Du_0);
        }
        
        
        STATE Du_n = (forceFrac_normal_Ef * Vm)/(forceFrac_normal_Ef+Kni*Vm);
        mem.SetFrac_normal(normal);
        mem.SetForceFrac_normal_n(forceFrac_normal_Ef);
        mem.SetDu_n(Du_n);
        mem.SetForceFrac_n(forceFrac);
        
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

template <class TMEM>
void TPZMatFractureBB<TMEM>::Correctnormal(TPZManVector<STATE,3>  &normal)
{
    REAL FracOrient = m_simulation_data->Get_FractureOrient().find(Frac_ID())->second;
    
    if(FracOrient==-1 && normal[1]<0.){
        normal[0]=-normal[0];
        normal[1]=-normal[1];
    }
    
    if(FracOrient==1 && normal[1]>0.){
        normal[0]=-normal[0];
        normal[1]=-normal[1];
    }
    int frac_Id = Frac_ID();

    if(frac_Id==10){
        normal[0]=-normal[0];
        normal[1]=-normal[1];
    }
    
    if(frac_Id==14){
        normal[0]=-normal[0];
        normal[1]=-normal[1];
    }
    
}

template class TPZMatFractureBB<TPZMemoryFracDFN>;

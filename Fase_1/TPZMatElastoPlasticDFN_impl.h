//
//  TPZMatElastoPlasticDFN_impl.h
//  Benchmark0a
//
//  Created by Omar Durán on 9/18/18.
//


#include "TPZMatElastoPlasticDFN.h"
#include "pzlog.h"
#include "pzbndcond.h"


#ifdef LOG4CXX
static LoggerPtr elastoplasticLogger(Logger::getLogger("pz.material.pzElastoPlasticDFN"));
static LoggerPtr updatelogger(Logger::getLogger("pz.material.pzElastoPlastic.update"));
static LoggerPtr ceckconvlogger(Logger::getLogger("checkconvmaterial"));
#endif

template <class T, class TMEM>
TPZMatElastoPlasticDFN<T,TMEM>::TPZMatElastoPlasticDFN() : TPZMatWithMem<TMEM>(), fForce(), fRhoB(0), fPostProcessDirection(), fTol(1.e-6)
{
    fForce.Resize(3,0);
    fForce[1] = 0.; // proper gravity acceleration in m/s^2
    fPostProcessDirection.Resize(3,0);
    fPostProcessDirection[0] = 1.;
    m_simulation_data = NULL;
    
#ifdef LOG4CXX
    if(elastoplasticLogger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>() constructor called ***";
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
#endif
    
}

template <class T, class TMEM>
TPZMatElastoPlasticDFN<T,TMEM>::TPZMatElastoPlasticDFN(int id, int dim) : TPZMatWithMem<TMEM>(id), fForce(), fRhoB(0), fPostProcessDirection(), fTol(1.e-6)
{
    fForce.Resize(3,0);
    fForce[1] = 0.; // proper gravity acceleration in m/s^2 -> 1=y 0=x 2=z
    fPostProcessDirection.Resize(3,0);
    fPostProcessDirection[0] = 1.;
    m_simulation_data = NULL;
    TPZPlasticState<STATE> def;
    
    if(dim!=2){
        DebugStop();
    }
    
    fDimension = dim;
    
#ifdef LOG4CXX
    if (elastoplasticLogger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>(int id) constructor called with id = " << id << " ***";
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
#endif
    
}

template <class T, class TMEM>
TPZMatElastoPlasticDFN<T,TMEM>::TPZMatElastoPlasticDFN(const TPZMatElastoPlasticDFN &mat) : TPZMatWithMem<TMEM>(mat),
fForce(mat.fForce), fRhoB(mat.fRhoB), fPostProcessDirection(mat.fPostProcessDirection),
fPlasticity(mat.fPlasticity), fTol(mat.fTol)
{
    m_simulation_data       = mat.m_simulation_data;
    
#ifdef LOG4CXX
    if(elastoplasticLogger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>() copy constructor called ***";
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
#endif
}


template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T,TMEM>::SetPlasticity(T & plasticity)
{
#ifdef LOG4CXX
    if(elastoplasticLogger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>::SetUpPlasticity ***";
        sout << "\n with plasticity argument:\n";
        plasticity.Print(sout);
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
#endif
    
    fPlasticity = plasticity;
    
    //fPlasticity.SetTensionSign(1);
    
    T plastloc(fPlasticity);
    
    TMEM memory;
    
    memory.GetPlasticState_n() = plastloc.GetState();
    
    plastloc.ApplyStrainComputeSigma(memory.GetPlasticState_n().fEpsT, memory.GetSigma());
    
    this->SetDefaultMem(memory);
    
#ifdef LOG4CXX
    if(elastoplasticLogger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "<< TPZMatElastoPlastic<T,TMEM>::SetUpPlasticity ***";
        sout << "\n with computed stresses:\n";
        sout << memory.GetSigma();
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
#endif
    
}


template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T,TMEM>::SetBulkDensity(REAL & RhoB)
{
    fRhoB = RhoB;
}

template <class T, class TMEM>
TPZMatElastoPlasticDFN<T,TMEM>::~TPZMatElastoPlasticDFN()
{
    
}

template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T,TMEM>::SetSimulationData(TPZSimulationData * simulation_data){
    m_simulation_data = simulation_data;
}


template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T,TMEM>::Print(std::ostream &out, const int memory)
{
    out << this->Name();
    out << "\n with template argurment T = " << fPlasticity.Name();
    out << "\n Base material Data:\n";
    TPZMatWithMem<TMEM>::PrintMem(out, memory);
    out << "\n Localy defined members:";
    out << "\n Body Forces: " << fForce;
    out << "\n Post process direction: " << fPostProcessDirection;
    out << "\n Tolerance for internal post processing iterations: " << fTol;
    out << "\n Internal plasticity <T> member:\n";
    fPlasticity.Print(out);
}

template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T,TMEM>::Print(std::ostream &out)
{
    out << this->Name();
    out << "\n with template argurment T = " << fPlasticity.Name();
    out << "\n Localy defined members:";
    out << "\n Body Forces: " << fForce;
    out << "\n Post process direction: " << fPostProcessDirection;
    out << "\n Tolerance for internal post processing iterations: " << fTol;
    out << "\n Internal plasticity <T> member:\n";
    fPlasticity.Print(out);
}

template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T,TMEM>::ApplyDeltaStrain(TPZMaterialData & data, TPZFMatrix<REAL> & DeltaStrain,TPZFMatrix<REAL> & Stress)
{
    
    int intPt = data.intGlobPtIndex;
    T plasticloc(fPlasticity);
    
    plasticloc.SetState(this->MemItem(intPt).GetPlasticState_n());
    
 //   UpdateMaterialCoeficients(data.x,plasticloc);
    TPZTensor<REAL> EpsT, Sigma;
    
    EpsT.CopyFrom(DeltaStrain);
    EpsT.Add(plasticloc.GetState().fEpsT, 1.);
    
    plasticloc.ApplyStrainComputeSigma(EpsT, Sigma);
    Sigma.CopyTo(Stress);
    
    if(TPZMatWithMem<TMEM>::fUpdateMem)
    {
        this->MemItem(intPt).GetSigma_n()        = Sigma;
        this->MemItem(intPt).GetPlasticState_n() = plasticloc.GetState();
      //  this->MemItem(intPt).fPlasticSteps = plasticloc.IntegrationSteps();
        int solsize = data.sol[0].size();
        for(int i=0; i<solsize; i++)
        {
            this->MemItem(intPt).Getu_n()[i] += data.sol[0][i];
        }
    }
    
}


template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T,TMEM>::ApplyDeltaStrainComputeDep(TPZMaterialData & data, TPZFMatrix<REAL> & DeltaStrain,TPZFMatrix<REAL> & Stress, TPZFMatrix<REAL> & Dep)
{

    int intPt = data.intGlobPtIndex;
    T plasticloc(fPlasticity);
    
    plasticloc.SetState(this->MemItem(intPt).GetPlasticState_n());
    
//    UpdateMaterialCoeficients(data.x,plasticloc);
    TPZTensor<REAL> EpsT, Sigma;
    
    EpsT.CopyFrom(DeltaStrain);
    EpsT.Add(plasticloc.GetState().fEpsT, 1.);
    
    plasticloc.ApplyStrainComputeSigma(EpsT, Sigma, &Dep);
    Sigma.CopyTo(Stress);
    
    if(TPZMatWithMem<TMEM>::fUpdateMem)
    {
        this->MemItem(intPt).GetSigma_n()        = Sigma;
        this->MemItem(intPt).GetPlasticState_n() = plasticloc.GetState();
 //       this->MemItem(intPt).fPlasticSteps = plasticloc.IntegrationSteps();
        int solsize = data.sol[0].size();
        for(int i=0; i<solsize; i++)
        {
            this->MemItem(intPt).Getu_n()[i] += data.sol[0][i];
        }
    }
    
    
}

template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T, TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef) {
    
    TPZFMatrix<REAL> &dphi = data.dphix, dphiXY;
    TPZFMatrix<REAL> &phi = data.phi;
    TPZFMatrix<REAL> &axes = data.axes, axesT;
    
    axes.Transpose(&axesT);
    axesT.Multiply(dphi, dphiXY);
    
    const int phr = phi.Rows();
    
    TPZFNMatrix<4> Deriv(2, 2);
    TPZFNMatrix<36> Dep(6, 6,0.0);
    TPZFNMatrix<6> DeltaStrain(6, 1);
    TPZFNMatrix<6> Stress(6, 1);
    int ptindex = data.intGlobPtIndex;
    
    if (TPZMatWithMem<TMEM>::fUpdateMem && data.sol.size() > 0) {
        // Loop over the solutions if update memory is true
        TPZSolVec locsol(data.sol);
        TPZGradSolVec locdsol(data.dsol);
        int numsol = locsol.size();
        
        for (int is = 0; is < numsol; is++) {
            data.sol[0] = locsol[is];
            data.dsol[0] = locdsol[is];
            
            this->ComputeDeltaStrainVector(data, DeltaStrain);
            this->ApplyDeltaStrainComputeDep(data, DeltaStrain, Stress, Dep);
        }
    } else {
        this->ComputeDeltaStrainVector(data, DeltaStrain);
        this->ApplyDeltaStrainComputeDep(data, DeltaStrain, Stress, Dep);
    }
    
#ifdef MACOS
    feclearexcept(FE_ALL_EXCEPT);
    if (fetestexcept(/*FE_DIVBYZERO*/ FE_ALL_EXCEPT)) {
        std::cout << "division by zero reported\n";
        DebugStop();
    }
#endif
    
#ifdef LOG4CXX
    if (elastoplasticLogger->isDebugEnabled()) {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>::Contribute ***";
        sout << "\nIntegration Local Point index = " << data.intGlobPtIndex;
        sout << "\nIntegration Global Point index = " << data.intGlobPtIndex;
        sout << "\ndata.axes = " << data.axes;
        sout << "\nDep " << endl;
        sout << Dep(_XX_, _XX_) << "\t" << Dep(_XX_, _YY_) << "\t" << Dep(_XX_, _XY_) << "\n";
        sout << Dep(_YY_, _XX_) << "\t" << Dep(_YY_, _YY_) << "\t" << Dep(_YY_, _XY_) << "\n";
        sout << Dep(_XY_, _XX_) << "\t" << Dep(_XY_, _YY_) << "\t" << Dep(_XY_, _XY_) << "\n";
        
        sout << "\nStress " << endl;
        sout << Stress(_XX_, 0) << "\t" << Stress(_YY_, 0) << "\t" << Stress(_XY_, 0) << "\n";
        
        sout << "\nDELTA STRAIN " << endl;
        sout << DeltaStrain(0, 0) << "\t" << DeltaStrain(1, 0) << "\t" << DeltaStrain(2, 0) << "\n";
        sout << "data.phi" << data.phi;
        
        LOGPZ_DEBUG(elastoplasticLogger, sout.str().c_str());
    }
#endif
    ptindex = 0;
    int nstate = NStateVariables();
    REAL val;
    
    TPZManVector<STATE, 5> ForceLoc(this->fForce);
    if (this->fForcingFunction) {
        this->fForcingFunction->Execute(data.x, ForceLoc);
    }
    
    int in;
    for (in = 0; in < phr; in++) {
        
        val = ForceLoc[0] * phi(in, 0);
        val += Stress(_XX_, 0) * dphiXY(0, in);
        val += Stress(_XY_, 0) * dphiXY(1, in);
        ef(in * nstate + 0, 0) += weight * val;
        
        val = ForceLoc[1] * phi(in, 0);
        val += Stress(_XY_, 0) * dphiXY(0, in);
        val += Stress(_YY_, 0) * dphiXY(1, in);
        ef(in * nstate + 1, 0) += weight * val;
        
        for (int jn = 0; jn < phr; jn++) {
            for (int ud = 0; ud < 2; ud++) {
                for (int vd = 0; vd < 2; vd++) {
                    Deriv(vd, ud) = dphiXY(vd, in) * dphiXY(ud, jn);
                }
            }
            
            val = 2. * Dep(_XX_, _XX_) * Deriv(0, 0); //dvdx*dudx
            val += Dep(_XX_, _XY_) * Deriv(0, 1); //dvdx*dudy
            val += 2. * Dep(_XY_, _XX_) * Deriv(1, 0); //dvdy*dudx
            val += Dep(_XY_, _XY_) * Deriv(1, 1); //dvdy*dudy
            val *= 0.5;
            ek(in * nstate + 0, jn * nstate + 0) += weight * val;
            
            val = Dep(_XX_, _XY_) * Deriv(0, 0);
            val += 2. * Dep(_XX_, _YY_) * Deriv(0, 1);
            val += Dep(_XY_, _XY_) * Deriv(1, 0);
            val += 2. * Dep(_XY_, _YY_) * Deriv(1, 1);
            val *= 0.5;
            ek(in * nstate + 0, jn * nstate + 1) += weight * val;
            
            val = 2. * Dep(_XY_, _XX_) * Deriv(0, 0);
            val += Dep(_XY_, _XY_) * Deriv(0, 1);
            val += 2. * Dep(_YY_, _XX_) * Deriv(1, 0);
            val += Dep(_YY_, _XY_) * Deriv(1, 1);
            val *= 0.5;
            ek(in * nstate + 1, jn * nstate + 0) += weight * val;
            
            val = Dep(_XY_, _XY_) * Deriv(0, 0);
            val += 2. * Dep(_XY_, _YY_) * Deriv(0, 1);
            val += Dep(_YY_, _XY_) * Deriv(1, 0);
            val += 2. * Dep(_YY_, _YY_) * Deriv(1, 1);
            val *= 0.5;
            ek(in * nstate + 1, jn * nstate + 1) += weight * val;
        }
    }
    
#ifdef LOG4CXX
    if (elastoplasticLogger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "<<< TPZMatElastoPlastic2D<T,TMEM>::Contribute ***";
        sout << " Resultant rhs vector:\n" << ef;
        sout << " Resultant stiff vector:\n" << ek;
        LOGPZ_DEBUG(elastoplasticLogger, sout.str().c_str());
    }
#endif
}

template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T, TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef) {
    TPZFMatrix<REAL> &dphi = data.dphix;
    TPZFMatrix<REAL> &phi = data.phi;
    TPZFMatrix<REAL> &axes = data.axes;
    TPZFNMatrix<9, REAL> axesT;
    TPZFNMatrix<50, REAL> dphiXY;
    
    axes.Transpose(&axesT);
    axesT.Multiply(dphi, dphiXY);
    
    const int phr = phi.Rows();
    
    //TPZFNMatrix<36> Deriv(6, 6);
    TPZFNMatrix<6> DeltaStrain(6, 1);
    TPZFNMatrix<6> Stress(6, 1);
    int ptindex = data.intGlobPtIndex;
    
    if (TPZMatWithMem<TMEM>::fUpdateMem && data.sol.size() > 0) {
        // Loop over the solutions if update memory is true
        //TPZFNMatrix<9> Dep(3, 3);
        
        TPZSolVec locsol(data.sol);
        TPZGradSolVec locdsol(data.dsol);
        int numsol = locsol.size();
        
        for (int is = 0; is < numsol; is++) {
            data.sol[0] = locsol[is];
            data.dsol[0] = locdsol[is];
            
            this->ComputeDeltaStrainVector(data, DeltaStrain);
            this->ApplyDeltaStrain(data, DeltaStrain, Stress);
            //this->ApplyDeltaStrainComputeDep(data, DeltaStrain, Stress, Dep);
        }
    } else {
        this->ComputeDeltaStrainVector(data, DeltaStrain);
        this->ApplyDeltaStrain(data, DeltaStrain, Stress);
        //        this->ApplyDeltaStrainComputeDep(data, DeltaStrain, Stress, Dep);
    }
#ifdef MACOS
    feclearexcept(FE_ALL_EXCEPT);
    if (fetestexcept(/*FE_DIVBYZERO*/ FE_ALL_EXCEPT)) {
        std::cout << "division by zero reported\n";
        DebugStop();
    }
#endif
    
#ifdef LOG4CXX
    if (elastoplasticLogger->isDebugEnabled()) {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>::Contribute ***";
        sout << "\nIntegration Local Point index = " << data.intGlobPtIndex;
        sout << "\nIntegration Global Point index = " << data.intGlobPtIndex;
        sout << "\ndata.axes = " << data.axes;
        sout << "\nStress " << endl;
        sout << Stress(_XX_, 0) << "\t" << Stress(_YY_, 0) << "\t" << Stress(_XY_, 0) << "\n";
        sout << "\nDELTA STRAIN " << endl;
        sout << DeltaStrain(0, 0) << "\t" << DeltaStrain(1, 0) << "\t" << DeltaStrain(2, 0) << "\n";
        sout << "data.phi" << data.phi;
        
        LOGPZ_DEBUG(elastoplasticLogger, sout.str().c_str());
    }
#endif
    ptindex = 0;
    int nstate = NStateVariables();
    REAL val;
    
    TPZManVector<STATE, 3> ForceLoc(this->fForce);
    if (this->fForcingFunction) {
        this->fForcingFunction->Execute(data.x, ForceLoc);
    }
    
    int in;
    for (in = 0; in < phr; in++) {
        val = ForceLoc[0] * phi(in, 0);
        val += Stress(_XX_, 0) * dphiXY(0, in);
        val += Stress(_XY_, 0) * dphiXY(1, in);
        ef(in * nstate + 0, 0) += weight * val;
        
        val = ForceLoc[1] * phi(in, 0);
        val += Stress(_XY_, 0) * dphiXY(0, in);
        val += Stress(_YY_, 0) * dphiXY(1, in);
        ef(in * nstate + 1, 0) += weight * val;
    }
    
    
#ifdef LOG4CXX
    if (elastoplasticLogger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "<<< TPZMatElastoPlastic2D<T,TMEM>::Contribute ***";
        sout << " Resultant rhs vector:\n" << ef;
        LOGPZ_DEBUG(elastoplasticLogger, sout.str().c_str());
    }
#endif
    
}

//template <class T, class TMEM>
//void TPZMatElastoPlasticDFN<T,TMEM>::FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
//{
//
//
//    //TPZMatWithMem<TMEM>::FillBoundaryConditionDataRequirement(type,data);
//    data.fNeedsSol = true;
//    if (type == 4 || type ==5 || type == 6) {
//        data.fNeedsNormal = true;
//    }
//    else {
//        data.fNeedsNormal = false;
//    }
//}

template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T, TMEM>::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc) {
  
    TPZBndCondWithMem<TPZMemoryBCDFN> & bc_with_memory = dynamic_cast<TPZBndCondWithMem<TPZMemoryBCDFN> &>(bc);
    int gp_index = data.intGlobPtIndex;
    if (m_simulation_data->Get_must_accept_solution_Q()) {
        
//        if (m_simulation_data->GetTransferCurrentToLastQ()) {
//            bc_with_memory.MemItem(gp_index).Setu(bc_with_memory.MemItem(gp_index).Getu_n()) ;
//            return;
//        }
        
        
            TPZManVector<STATE,3> delta_u    = data.sol[0];
            TPZManVector<STATE,3> u_n(fDimension,0.0);
            TPZManVector<STATE,3> u(bc_with_memory.MemItem(gp_index).Getu());
            for (int i = 0; i < fDimension; i++) {
                u_n[i] = delta_u[i] + u[i];
            }
            bc_with_memory.MemItem(gp_index).Setu_n(u_n);
            
    }
    
    TPZFMatrix<REAL>  &phi = data.phi;
    TPZManVector<STATE,3> delta_u    = data.sol[0];
    TPZManVector<STATE,3> u_n(fDimension,0.0);
    TPZManVector<STATE,3> u(bc_with_memory.MemItem(gp_index).Getu_n());
    for (int i = 0; i < fDimension; i++) {
        u_n[i] = delta_u[i] + u[i];
    }
    
    
    int phr = phi.Rows();
    int in, jn, idf, jdf;
    REAL BigNumber = TPZDiscontinuousGalerkin::gBigNumber;
    int nstate = NStateVariables();

    REAL v2[2];
    v2[0] = bc.Val2()(0, 0);
    v2[1] = bc.Val2()(1, 0);
    
    
    TPZFMatrix<REAL> &v1 = bc.Val1();
    switch (bc.Type()) {
        case 0: // Dirichlet condition
            for (in = 0; in < phr; in++) {
                ef(nstate * in + 0, 0) += BigNumber * (data.sol[0][0]-v2[0]) * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += BigNumber * (data.sol[0][1]-v2[1]) * phi(in, 0) * weight;
                
                for (jn = 0; jn < phr; jn++) {
                    ek(nstate * in + 0, nstate * jn + 0) += BigNumber * phi(in, 0) * phi(jn, 0) * weight;
                    ek(nstate * in + 1, nstate * jn + 1) += BigNumber * phi(in, 0) * phi(jn, 0) * weight;
                    
                }//jn
            }//in
            break;
            
        case 1: // Neumann condition
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) -= v2[0] * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) -= v2[1] * phi(in, 0) * weight;
            }
            break;
            
        case 2: // Mixed condition
        {
            TPZFNMatrix<2, STATE> res(2, 1, 0.);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    res(i, 0) += bc.Val1()(i, j) * data.sol[0][j];
                }
            }
            
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) += (v2[0] - res(0, 0)) * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += (v2[1] - res(1, 0)) * phi(in, 0) * weight;
                for (jn = 0; jn < phi.Rows(); jn++) {
                    for (idf = 0; idf < 2; idf++) {
                        for (jdf = 0; jdf < 2; jdf++) {
                            ek(nstate * in + idf, nstate * jn + jdf) += bc.Val1()(idf, jdf) * phi(in, 0) * phi(jn, 0) * weight;
                            //BUG FALTA COLOCAR VAL2
                            //DebugStop();
                        }
                    }
                }
            }//in
        }
            break;
            
        case 3: // Directional Null Dirichlet - displacement is set to null in the non-null vector component direction
            for (in = 0; in < phr; in++) {
                ef(nstate * in + 0, 0) += BigNumber * (u_n[0] - 0.0) * v2[0] * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += BigNumber * (u_n[1] - 0.0) * v2[1] * phi(in, 0) * weight;
                for (jn = 0; jn < phr; jn++) {
                    ek(nstate * in + 0, nstate * jn + 0) += BigNumber * phi(in, 0) * phi(jn, 0) * weight * v2[0];
                    ek(nstate * in + 1, nstate * jn + 1) += BigNumber * phi(in, 0) * phi(jn, 0) * weight * v2[1];
                }//jn
            }//in
            break;
        
            
        case 4: // stressField Neumann condition
            v2[0] = v1(0, 0) * data.normal[0] + v1(0, 1) * data.normal[1];
            v2[1] = v1(1, 0) * data.normal[0] + v1(1, 1) * data.normal[1];
            // The normal vector points towards the neighbor. The negative sign is there to
            // reflect the outward normal vector.
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) += v2[0] * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += v2[1] * phi(in, 0) * weight;
            }
            break;
        case 5://PRESSAO DEVE SER POSTA NA POSICAO 0 DO VETOR v2
        {
            TPZFNMatrix<2, STATE> res(2, 1, 0.);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    res(i, 0) += data.normal[i] * bc.Val1()(i, j) * data.sol[0][j] * data.normal[j];
                }
            }
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) += (v2[0] * data.normal[0] - res(0, 0)) * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += (v2[0] * data.normal[1] - res(1, 0)) * phi(in, 0) * weight;
                for (jn = 0; jn < phi.Rows(); jn++) {
                    for (idf = 0; idf < 2; idf++) {
                        for (jdf = 0; jdf < 2; jdf++) {
                            ek(nstate * in + idf, nstate * jn + jdf) += bc.Val1()(idf, jdf) * data.normal[idf] * data.normal[jdf] * phi(in, 0) * phi(jn, 0) * weight;
                            // BUG FALTA COLOCAR VAL2
                            // DebugStop();
                        }
                    }
                }
                
            }
        }
            break;
            
        case 6://PRESSAO DEVE SER POSTA NA POSICAO 0 DO VETOR v2
        {
            TPZFNMatrix<2, STATE> res(2, 1, 0.);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    res(i, 0) += bc.Val1()(i, j) * data.sol[0][j];
                }
            }
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) += (v2[0] * data.normal[0] - res(0, 0)) * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += (v2[0] * data.normal[1] - res(1, 0)) * phi(in, 0) * weight;
                for (jn = 0; jn < phi.Rows(); jn++) {
                    for (idf = 0; idf < 2; idf++) {
                        for (jdf = 0; jdf < 2; jdf++) {
                            ek(nstate * in + idf, nstate * jn + jdf) += bc.Val1()(idf, jdf) * phi(in, 0) * phi(jn, 0) * weight;
                            // BUG FALTA COLOCAR VAL2
                            // DebugStop();
                        }
                    }
                }
                
            }
            
        }
            break;
            
        default:
#ifdef LOG4CXX
            if (elastoplasticLogger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "<<< TPZMatElastoPlastic2D<T,TMEM>::ContributeBC *** WRONG BOUNDARY CONDITION TYPE = " << bc.Type();
                LOGPZ_ERROR(elastoplasticLogger, sout.str().c_str());
            }
#endif
            PZError << "TPZMatElastoPlastic2D::ContributeBC error - Wrong boundary condition type" << std::endl;
    }
}

template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T, TMEM>::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
    TPZFMatrix<REAL> ek_fake;
    ek_fake.Resize(ef.Rows(),ef.Rows());
    this->ContributeBC(data, weight, ek_fake, ef, bc);
}



template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T,TMEM>::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
{
 
    long gp_index = data.intGlobPtIndex;
    TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    Solout.Resize( this->NSolutionVariables(var));
    
    TPZTensor<REAL> epsilon_t = memory.GetPlasticState_n().fEpsT;
    TPZTensor<REAL> epsilon_p = memory.GetPlasticState_n().fEpsP;
    
    switch (var) {
        case 0:
        {
            Solout[0] = memory.Getu_n()[0];
        }
            break;
        case 1:
        {
            Solout[0] = memory.Getu_n()[1];
        }
            break;
        case 2:
        {
            Solout[0] = memory.Getu_n()[2];
        }
            break;
        case 3:
        {
            Solout[0] = memory.GetSigma_n().XX();
        }
            break;
        case 4:
        {
            Solout[0] = memory.GetSigma_n().XY();
        }
            break;
        case 5:
        {
            Solout[0] = memory.GetSigma_n().XZ();
        }
            break;
        case 6:
        {
            Solout[0] = memory.GetSigma_n().YY();
        }
            break;
        case 7:
        {
            Solout[0] = memory.GetSigma_n().YZ();
        }
            break;
        case 8:
        {
            Solout[0] = memory.GetSigma_n().ZZ();
        }
            break;
        case 9:
        {
            Solout[0] = epsilon_t.XX();
        }
            break;
        case 10:
        {
            Solout[0] = epsilon_t.XY();
        }
            break;
        case 11:
        {
            Solout[0] = epsilon_t.XZ();
        }
            break;
        case 12:
        {
            Solout[0] = epsilon_t.YY();
        }
            break;
        case 13:
        {
            Solout[0] = epsilon_t.YZ();
        }
            break;
        case 14:
        {
            Solout[0] = epsilon_t.ZZ();
        }
            break;
        case 15:
        {
            Solout[0] = epsilon_p.XX();
        }
            break;
        case 16:
        {
            Solout[0] = epsilon_p.XY();
        }
            break;
        case 17:
        {
            Solout[0] = epsilon_p.XZ();
        }
            break;
        case 18:
        {
            Solout[0] = epsilon_p.YY();
        }
            break;
        case 19:
        {
            Solout[0] = epsilon_p.YZ();
        }
            break;
        case 20:
        {
            Solout[0] = epsilon_p.ZZ();
        }
            break;
        default:
        {
            std::cout << "TPMRSElastoPlastic<T,TMEM>:: Variable not implemented." << std::endl;
            DebugStop();
        }
            break;
    }
 
    
}

template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T, TMEM>::ComputeDeltaStrainVector(TPZMaterialData & data, TPZFMatrix<REAL> &DeltaStrain) {
    TPZFNMatrix<9> DSolXYZ(3, 3, 0.);
    data.axes.Multiply(data.dsol[0], DSolXYZ, 1/*transpose*/);
    if (DeltaStrain.Rows() != 6) {
        DebugStop();
    }
    //  DeltaStrain.Redim(3,1);
    DeltaStrain(_XX_, 0) = DSolXYZ(0, 0);
    DeltaStrain(_YY_, 0) = DSolXYZ(1, 1);
    DeltaStrain(_XY_, 0) = 0.5 * (DSolXYZ(1, 0) + DSolXYZ(0, 1));
    DeltaStrain(_XZ_, 0) = 0.;
    DeltaStrain(_YZ_, 0) = 0.;
    DeltaStrain(_ZZ_, 0) = 0.;
}

template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T,TMEM>::EigenValues(TPZFMatrix<REAL> & vectorTensor, TPZVec<REAL> & ev)
{
    TPZFNMatrix<9> Tensor(3,3);
    ev.Resize(3);
    this->vectorToTensor(vectorTensor, Tensor);
    int64_t numiterations = 1000;
    
#ifdef PZDEBUG
    bool result = Tensor.SolveEigenvaluesJacobi(numiterations, fTol, &ev);
    if (result == false){
        PZError << __PRETTY_FUNCTION__ << " - ERROR! - result = false - numiterations = " << numiterations << " - tol = " << fTol << std::endl;
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "<<< TPZMatElastoPlasticDFN<T,TMEM>::EigenValues *** not solved within " << numiterations << " iterations";
            sout << "\n vectorTensor = " << vectorTensor;
            LOGPZ_ERROR(elastoplasticLogger,sout.str().c_str());
        }
#endif
    }
#else
    Tensor.SolveEigenvaluesJacobi(numiterations, fTol, &ev);
#endif
}

template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T,TMEM>::EigenVectors(TPZFMatrix<REAL> &vectorTensor, TPZVec< REAL > &Solout, int direction)
{
    TPZFNMatrix<9> Tensor(3,3);
    this->vectorToTensor(vectorTensor, Tensor);
    
    TPZManVector<REAL,3> Eigenvalues(3);
    TPZFNMatrix<9> Eigenvectors(3,3);
    
    int64_t numiterations = 1000;
#ifdef PZDEBUG
    bool result = Tensor.SolveEigensystemJacobi(numiterations, fTol, Eigenvalues, Eigenvectors);
    if (result == false){
        PZError << __PRETTY_FUNCTION__ << " - ERROR! - result = false - numiterations = " << numiterations << " - tol = " << fTol << std::endl;
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "<<< TPZMatElastoPlasticDFN<T,TMEM>::EigenVectors *** not solved within " << numiterations << " iterations";
            sout << "\n vectorTensor = " << vectorTensor;
            LOGPZ_ERROR(elastoplasticLogger,sout.str().c_str());
        }
#endif
    }
#else
    Tensor.SolveEigensystemJacobi(numiterations, fTol, Eigenvalues, Eigenvectors);
#endif
    Solout.Resize(3);
    for(int i = 0; i < 3; i++) Solout[i] = Eigenvectors(direction,i);
}

template <class T, class TMEM>
TPZMaterial * TPZMatElastoPlasticDFN<T,TMEM>::NewMaterial()
{
    return new TPZMatElastoPlasticDFN<T,TMEM>(*this);
}
/*
 void TPZMatElastoPlasticDFN::SetData(std::istream &data)
 {
 TPZMaterial::SetData(data);
 data >> fDeltaT; // to be removed in the elastoplastic material and readded to the poroelastoplastic material
 }*/

#include "TPZSandlerExtended.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"

template <class T, class TMEM>
std::string TPZMatElastoPlasticDFN<T,TMEM>::Name() {
    return "TPZMatElastoPlasticDFN<T,TMEM>";
}

template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T, TMEM>::Write(TPZStream &buf, int withclassid) const {
    TPZMatWithMem<TMEM>::Write(buf, withclassid);
    
    buf.Write(&fForce[0], 3);
    buf.Write(&fPostProcessDirection[0], 3);
    fPlasticity.Write(buf, withclassid);
    buf.Write(&fTol, 1);
}

template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T, TMEM>::Read(TPZStream &buf, void *context) {
    //    TPZSavable::Read(buf, context);
    
    TPZMatWithMem<TMEM>::Read(buf, context);
    
    buf.Read(&fForce[0], 3);
    buf.Read(&fPostProcessDirection[0], 3);
    fPlasticity.Read(buf, context);
    buf.Read(&fTol, 1);
}

template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T,TMEM>::SetTol(const REAL & tol)
{
    fTol = tol;
}

template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T,TMEM>::SetBulkDensity(const REAL & bulk)
{
    fRhoB = bulk;
}

template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T,TMEM>::vectorToTensor(const TPZFMatrix<REAL> & vectorTensor, TPZFMatrix<REAL> & Tensor)
{
    TPZTensor<REAL> vecT;
    vecT.CopyFrom(vectorTensor);
    vecT.CopyToTensor(Tensor);
}

template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T,TMEM>::FillDataRequirements(TPZMaterialData &data){
    
    TPZMatWithMem<TMEM>::FillDataRequirements(data);
    
    data.fNeedsSol = true;
    data.fNeedsNormal = false;
    data.fNeedsHSize = false;
    data.fNeedsNeighborCenter = false;
}

template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T,TMEM>::FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
{
    
    
    //TPZMatWithMem<TMEM>::FillBoundaryConditionDataRequirement(type,data);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
}

template <class T, class TMEM>
void TPZMatElastoPlasticDFN<T,TMEM>::Errors(TPZVec<REAL> &x,TPZVec<REAL> &u, TPZFMatrix<REAL> &dudx,
                                         TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux,
                                         TPZVec<REAL> &u_exact,TPZFMatrix<REAL> &du_exact,TPZVec<REAL> &values)
{
    int i, j;
    
    /** L2 norm */
    REAL L2 = 0.;
    for(i = 0; i < 3; i++) L2 += (u[i] - u_exact[i]) * (u[i] - u_exact[i]);
    
    /** H1 semi-norm */
    REAL SemiH1 = 0.;
    for(i = 0; i < 3; i++) for(j = 0; j < 3; j++) SemiH1 += (dudx(i,j) - du_exact(i,j)) * (dudx(i,j) - du_exact(i,j));
    
    /** H1 norm */
    REAL H1 = L2 + SemiH1;
    
    //values[1] : eror em norma L2
    values[1]  = L2;
    
    //values[2] : erro em semi norma H1
    values[2] = SemiH1;
    
    //values[0] : erro em norma H1 <=> norma Energia
    values[0]  = H1;
    
}

template <class T, class TMEM>
int TPZMatElastoPlasticDFN<T,TMEM>::VariableIndex(const std::string &name)
{
    if (!strcmp("ux", name.c_str())) return 0;
    if (!strcmp("uy", name.c_str())) return 1;
    if (!strcmp("uz", name.c_str())) return 2;
    if (!strcmp("sxx", name.c_str())) return 3;
    if (!strcmp("sxy", name.c_str())) return 4;
    if (!strcmp("sxz", name.c_str())) return 5;
    if (!strcmp("syy", name.c_str())) return 6;
    if (!strcmp("syz", name.c_str())) return 7;
    if (!strcmp("szz", name.c_str())) return 8;
    if (!strcmp("exx", name.c_str())) return 9;
    if (!strcmp("exy", name.c_str())) return 10;
    if (!strcmp("exz", name.c_str())) return 11;
    if (!strcmp("eyy", name.c_str())) return 12;
    if (!strcmp("eyz", name.c_str())) return 13;
    if (!strcmp("ezz", name.c_str())) return 14;
    if (!strcmp("epxx", name.c_str())) return 15;
    if (!strcmp("epxy", name.c_str())) return 16;
    if (!strcmp("epxz", name.c_str())) return 17;
    if (!strcmp("epyy", name.c_str())) return 18;
    if (!strcmp("epyz", name.c_str())) return 19;
    if (!strcmp("epzz", name.c_str())) return 20;
    return TPZMatWithMem<TMEM>::VariableIndex(name);
}

template <class T, class TMEM>
int TPZMatElastoPlasticDFN<T,TMEM>::NSolutionVariables(int var)
{
    switch(var) {
        case 0:
            return 1; // Scalar
        case 1:
            return 1; // Scalar
        case 2:
            return 1; // Scalar
        case 3:
            return 1; // Scalar
        case 4:
            return 1; // Scalar
        case 5:
            return 1; // Scalar
        case 6:
            return 1; // Scalar
        case 7:
            return 1; // Scalar
        case 8:
            return 1; // Scalar
        case 9:
            return 1; // Scalar
        case 10:
            return 1; // Scalar
        case 11:
            return 1; // Scalar
        case 12:
            return 1; // Scalar
        case 13:
            return 1; // Scalar
        case 14:
            return 1; // Scalar
        case 15:
            return 1; // Scalar
        case 16:
            return 1; // Scalar
        case 17:
            return 1; // Scalar
        case 18:
            return 1; // Scalar
        case 19:
            return 1; // Scalar
        case 20:
            return 1; // Scalar
    }
    return TPZMatWithMem<TMEM>::NSolutionVariables(var);
}



#include "meshgen.h"

#include "pzgmesh.h"
#include "pzgengrid.h"
#include "pzgeoel.h"
#include "TPZRefPatternTools.h"
#include "pzcheckgeom.h"
#include "TPZVTKGeoMesh.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"

#include "TPZMaterial.h"
#include "pzbuildmultiphysicsmesh.h"

#include "TPZSSpStructMatrix.h"

#ifndef USING_MKL
#include "pzskylstrmatrix.h"
#endif


#include "pzlog.h"

#ifdef _AUTODIFF
#include "fadType.h"

static FADFADREAL FADsin(FADFADREAL x)
{
    FADREAL_ sinaval = sin(x.val());
    FADREAL_ cosaval = cos(x.val());
    FADFADREAL sina(2,sinaval);
    for (int i=0; i<2; i++) {
        sina.fastAccessDx(i) = cosaval*x.dx(i);
    }
    return sina;
}

static FADFADREAL FADcos(FADFADREAL x)
{
    FADREAL_ sinaval = sin(x.val());
    FADREAL_ cosaval = cos(x.val());
    FADFADREAL cosa(2,cosaval);
    for (int i=0; i<2; i++) {
        cosa.fastAccessDx(i) = -sinaval*x.dx(i);
    }
    return cosa;
}

static FADFADREAL FADexp(FADFADREAL x)
{
    FADREAL_ expaval = exp(x.val());
    FADFADREAL expa(2,expaval);
    for (int i=0; i<2; i++) {
        expa.fastAccessDx(i) = expaval*x.dx(i);
    }
    return expa;
}

static FADFADREAL FADsqrt(FADFADREAL x)
{
    FADREAL_ fadres = sqrt(x.val());
    FADFADREAL resa(2,fadres);
    for (int i=0; i<2; i++) {
        resa.fastAccessDx(i) = REAL(0.5)/fadres*x.dx(i);
    }
    return resa;
}

static FADFADREAL FADatan(FADFADREAL x)
{
    FADREAL_ fadres = atan(x.val());
    FADFADREAL resa(2,fadres);
    for (int i=0; i<2; i++) {
        resa.fastAccessDx(i) = 1./(1+x.val()*x.val())*x.dx(i);
    }
    return resa;
}

template<class TVar>
void TElasticityExample1::uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp)
{
    disp[0] = TVar(1./27.)*x[0]*x[0]*x[1]*x[1]*cos(TVar(6.*M_PI)*x[0])*sin(TVar(7.*M_PI)*x[1]);

    disp[1] = TVar(0.2)*exp(x[1])*sin(TVar(4.*M_PI)*x[0]);
}

template<>
void TElasticityExample1::uxy(const TPZVec<FADFADREAL > &x, TPZVec<FADFADREAL > &disp)
{
    FADFADREAL tmp = (FADFADREAL)(1./27.)*x[0]*x[0]*x[1]*x[1];
    disp[0] = tmp*FADcos((FADFADREAL)(6.*M_PI)*x[0])*FADsin((FADFADREAL)(7.*M_PI)*x[1]);
    disp[1] = (FADFADREAL)(0.2)*FADexp(x[1])*FADsin((FADFADREAL)(4.*M_PI)*x[0]);
}

template<class TVar>
void TElasticityExample1::Elastic(const TPZVec<TVar> &x, TVar &Elast, TVar &nu)
{
    Elast = (TVar(100.) * (TVar(1.) + TVar(0.3) * sin(TVar(10 * M_PI) * (x[0] - TVar(0.5))) * cos(TVar(10. * M_PI) * x[1])));
//    Elast.val() = 1000.;
    nu = TVar(0.3);
}

template<>
void TElasticityExample1::Elastic(const TPZVec<double> &x, double &Elast, double &nu)
{
//    Elast = 1000.;
    Elast = (100. * (1. + 0.3 * sin(10 * M_PI * (x[0] - 0.5)) * cos(10. * M_PI * x[1])));
    nu = 0.3;
}

void TElasticityExample1::ElasticDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)
{
    TPZManVector<STATE> xstate(x.size());
    for (int i=0; i<xstate.size(); i++) {
        xstate[i] = x[i];
    }
    STATE E,nu;
    Elastic(xstate,E,nu);
    result[0] = E;
    result[1] = nu;
}


TPZAutoPointer<TPZFunction<STATE> > TElasticityExample1::ConstitutiveLawFunction()
{
    TPZAutoPointer<TPZFunction<STATE> > result;
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(TElasticityExample1::ElasticDummy,4);
    dummy->SetPolynomialOrder(4);
    result = TPZAutoPointer<TPZFunction<STATE> >(dummy);
    return result;
}

template<class TVar>
void TElasticityExample1::graduxy(const TPZVec<TVar> &x, TPZFMatrix<TVar> &grad)
{
    TPZManVector<Fad<TVar>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        Fad<TVar> temp = Fad<TVar>(2,i,x[i]);
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<TVar>,3> result(2);
    uxy(xfad,result);
    grad.Resize(2,2);
    for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++)
        {
            grad(i,j) = result[i].d(j);
        }
    }
//    for(int i=0; i<2; i++)
//    {
//        std::cout << "result " << result[i] << " dx " << result[i].dx() << std::endl;
//        for (int j=0; j<2; j++) {
//            std::cout << "i = " << i << " j = " << j << " grad " << grad(i,j) << " dx " << grad(i,j).dx() << std::endl;
//        }
//    }
}

void TElasticityExample1::GradU(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu)
{
    TPZManVector<Fad<REAL>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        Fad<REAL> temp = Fad<REAL>(2,i,x[i]);
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<REAL>,3> result(2);
    uxy(xfad,result);
    gradu.Redim(2,2);
    u[0] = result[0].val();
    u[1] = result[1].val();
    for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++)
        {
            gradu(i,j) = result[j].d(i);
        }
    }
    
}


template<>
void TElasticityExample1::graduxy(const TPZVec<Fad<REAL> > &x, TPZFMatrix<Fad<REAL> > &grad)
{
    TPZManVector<Fad<Fad<REAL> >,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        Fad<Fad<REAL> > temp = Fad<Fad<REAL> >(2,Fad<REAL>(2,0.));
        temp.val()= x[i];
        Fad<REAL> temp2(2,1.);
        temp.fastAccessDx(i) = temp2;
        xfad[i] = temp;
    }
    TPZManVector<Fad<Fad<REAL> >,3> result(2);
    uxy(xfad,result);
    grad.Resize(2,2);
    for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++)
        {
            grad(i,j) = result[i].d(j);
        }
    }
//    for(int i=0; i<2; i++)
//    {
//        std::cout << "result " << result[i] << " dx " << result[i].dx() << std::endl;
//        for (int j=0; j<2; j++) {
//            std::cout << "i = " << i << " j = " << j << " grad " << grad(i,j) << " dx " << grad(i,j).dx() << std::endl;
//        }
//    }
}

template<class TVar>
void TElasticityExample1::Sigma(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma)
{
    TPZFNMatrix<4,TVar> grad;
    TVar E, nu;
    Elastic(x, E, nu);
    TVar Fac = E/((TVar)1.+nu)/((TVar(1.)-TVar(2.)*nu));
    graduxy(x,grad);
    sigma.Resize(2,2);
    sigma(0,0) = Fac*((TVar(1.)-nu)*grad(0,0)+nu*grad(1,1));
    sigma(1,1) = Fac*((TVar(1.)-nu)*grad(1,1)+nu*grad(0,0));
    sigma(0,1) = E/(TVar(2.)*(TVar(1.)+nu))*(grad(0,1)+grad(1,0));
    sigma(1,0) = sigma(0,1);
//    for(int i=0; i<2; i++)
//    {
//        for (int j=0; j<2; j++) {
//            std::cout << "i = " << i << " j = " << j << " grad " << grad(i,j) << " sigma " << sigma(i,j) << std::endl;
//        }
//    }
}

template<>
void TElasticityExample1::Sigma(const TPZVec<Fad<REAL> > &x, TPZFMatrix<Fad<REAL> > &sigma)
{
    TPZFNMatrix<4,Fad<REAL> > grad;
    Fad<REAL>  E, nu;
    Elastic(x, E, nu);
    Fad<REAL>  Fac = E/(Fad<REAL>(1.)+nu)/((Fad<REAL>(1.)-Fad<REAL>(2.)*nu));
    graduxy(x,grad);
    sigma.Resize(2,2);
    sigma(0,0) = Fac*((Fad<REAL>(1.)-nu)*grad(0,0)+nu*grad(1,1));
    sigma(1,1) = Fac*((Fad<REAL>(1.)-nu)*grad(1,1)+nu*grad(0,0));
    sigma(0,1) = E/(Fad<REAL>(2.)*(Fad<REAL>(1.)+nu))*(grad(0,1)+grad(1,0));
    sigma(1,0) = sigma(0,1);
//    for(int i=0; i<2; i++)
//    {
//        for (int j=0; j<2; j++) {
//            std::cout << "i = " << i << " j = " << j << " grad " << grad(i,j) << " sigma " << sigma(i,j) << std::endl;
//        }
//    }
}

template<class TVar>
void TElasticityExample1::DivSigma(const TPZVec<TVar> &x, TPZVec<TVar> &divsigma)
{
    TPZManVector<Fad<TVar>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        xfad[i] = Fad<TVar>(2,i,x[i]);
    }
    TPZFNMatrix<4, Fad<TVar> > sigma(2,2);
    Sigma(xfad,sigma);
//    for(int i=0; i<2; i++)
//    {
//        for (int j=0; j<2; j++) {
//            std::cout << "i = " << i << " j = " << j <<  " sigma " << sigma(i,j) <<  " " << sigma(i,j).dx() << std::endl;
//        }
//    }

    divsigma[0] = sigma(0,0).dx(0)+sigma(0,1).dx(1);
    divsigma[1] = sigma(1,0).dx(0)+sigma(1,1).dx(1);
}

TPZAutoPointer<TPZFunction<STATE> > TElasticityExample1::ForcingFunction()
{
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Force,5);
    dummy->SetPolynomialOrder(5);
    TPZAutoPointer<TPZFunction<STATE> > result(dummy);
    return result;
}

TPZAutoPointer<TPZFunction<STATE> > TElasticityExample1::ValueFunction()
{
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(GradU,5);
    dummy->SetPolynomialOrder(5);
    TPZAutoPointer<TPZFunction<STATE> > result(dummy);
    return result;
    
}


template
void TElasticityExample1::DivSigma<REAL>(const TPZVec<REAL> &x, TPZVec<REAL> &divsigma);
template
void TElasticityExample1::Sigma<Fad<REAL> >(const TPZVec<Fad<REAL> > &x, TPZFMatrix<Fad<REAL> > &sigma);



template<class TVar>
void TLaplaceExample1::uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp)
{
    disp[0] = sin((TVar)M_PI*x[0])*sin((TVar)M_PI*x[1]);
    TVar r = sqrt(x[0]*x[0]+x[1]*x[1]);
    TVar atanco = (r-(TVar)0.5)*100.;
    TVar freq = 10.;
    TVar mult = (TVar(1)+TVar(0.3)*sin(TVar(M_PI)*x[0]*freq))*(TVar(1)+TVar(0.5)*cos(TVar(M_PI)*r*freq));
    disp[0] = atan(atanco)*mult;
}

template<>
void TLaplaceExample1::uxy(const TPZVec<FADFADREAL > &x, TPZVec<FADFADREAL > &disp)
{
//    disp[0] = FADsin((FADFADREAL)(M_PI)*x[0])*FADsin((FADFADREAL)(M_PI)*x[1]);
    FADFADREAL r = FADsqrt(x[0]*x[0]+x[1]*x[1]);
    FADFADREAL atanco = (r-(FADFADREAL)0.5)*100.;
    FADFADREAL freq = (FADFADREAL)10.;
    FADFADREAL mult = ((FADFADREAL)1.+(FADFADREAL)0.3*FADsin((FADFADREAL)M_PI*x[0]*freq))*((FADFADREAL)1.+(FADFADREAL)0.5*FADcos((FADFADREAL)M_PI*r*freq));
    disp[0] = FADatan(atanco)*mult;

}

template<class TVar>
void TLaplaceExample1::Permeability(const TPZVec<TVar> &x, TVar &Perm)
{
    Perm = (TVar)(1.);
}

void TLaplaceExample1::PermeabilityDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)
{
    TPZManVector<STATE,3> xloc(x.size());
    for (auto i : xloc) {
        xloc[i] = x[i];
    }
    STATE Perm;
    Permeability(xloc, Perm);
    deriv.Zero();
    deriv(0,0) = Perm;
    deriv(1,1) = Perm;
    deriv(2,0) = 1./Perm;
    deriv(3,1) = 1./Perm;
}

template<class TVar>
void TLaplaceExample1::graduxy(const TPZVec<TVar> &x, TPZVec<TVar> &grad)
{
    TPZManVector<Fad<TVar>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        Fad<TVar> temp = Fad<TVar>(2,i,x[i]);
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<TVar>,3> result(1);
    uxy(xfad,result);
    grad.resize(2);
    for (int i=0; i<2; i++)
    {
            grad[i] = result[0].d(i);
    }
    //    for(int i=0; i<2; i++)
    //    {
    //        std::cout << "result " << result[i] << " dx " << result[i].dx() << std::endl;
    //        for (int j=0; j<2; j++) {
    //            std::cout << "i = " << i << " j = " << j << " grad " << grad(i,j) << " dx " << grad(i,j).dx() << std::endl;
    //        }
    //    }
}

void TLaplaceExample1::GradU(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu)
{
    TPZManVector<Fad<REAL>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        Fad<REAL> temp = Fad<REAL>(2,i,x[i]);
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<REAL>,3> result(2);
    uxy(xfad,result);
    gradu.Redim(2,1);
    u[0] = result[0].val();
    for (int i=0; i<2; i++) {
        for (int j=0; j<1; j++)
        {
            gradu(i,j) = result[j].d(i);
        }
    }
    
}

template<class TVar>
void TLaplaceExample1::Sigma(const TPZVec<TVar> &x, TPZVec<TVar> &sigma)
{
    TPZManVector<TVar,3> grad;
    TVar Perm;
    Permeability(x, Perm);
    graduxy(x,grad);
    sigma.resize(2);
    sigma[0] = -Perm*grad[0];
    sigma[1] = -Perm*grad[1];
    
}

template<class TVar>
void TLaplaceExample1::DivSigma(const TPZVec<TVar> &x, TVar &divsigma)
{
    TPZManVector<Fad<TVar>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        xfad[i] = Fad<TVar>(2,i,x[i]);
    }
    TPZManVector<Fad<TVar>, 3> sigma(2);
    Sigma(xfad,sigma);
    //    for(int i=0; i<2; i++)
    //    {
    //        for (int j=0; j<2; j++) {
    //            std::cout << "i = " << i << " j = " << j <<  " sigma " << sigma(i,j) <<  " " << sigma(i,j).dx() << std::endl;
    //        }
    //    }
    
    divsigma = sigma[0].dx(0)+sigma[1].dx(1);
    
}


TPZAutoPointer<TPZFunction<STATE> > TLaplaceExample1::ForcingFunction()
{
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Force,5);
    dummy->SetPolynomialOrder(5);
    TPZAutoPointer<TPZFunction<STATE> > result(dummy);
    return result;
}

TPZAutoPointer<TPZFunction<STATE> > TLaplaceExample1::ValueFunction()
{
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(GradU,5);
    dummy->SetPolynomialOrder(5);
    TPZAutoPointer<TPZFunction<STATE> > result(dummy);
    return result;
    
}

TPZAutoPointer<TPZFunction<STATE> > TLaplaceExample1::ConstitutiveLawFunction()
{
    TPZAutoPointer<TPZFunction<STATE> > result;
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(PermeabilityDummy,4);
    dummy->SetPolynomialOrder(4);
    result = TPZAutoPointer<TPZFunction<STATE> >(dummy);
    return result;
    
}

//ExactFunc *Exact();


template<class TVar>
void TLaplaceExampleSmooth::uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp)
{
    disp[0] = x[0];
}

template<>
void TLaplaceExampleSmooth::uxy(const TPZVec<FADFADREAL > &x, TPZVec<FADFADREAL > &disp)
{
    disp[0] = x[0];
    
}

template<class TVar>
void TLaplaceExampleSmooth::Permeability(const TPZVec<TVar> &x, TVar &Perm)
{
    Perm = (TVar)(1.);
}

void TLaplaceExampleSmooth::PermeabilityDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)
{
    TPZManVector<STATE,3> xloc(x.size());
    for (auto i : xloc) {
        xloc[i] = x[i];
    }
    STATE Perm;
    Permeability(xloc, Perm);
    deriv.Zero();
    deriv(0,0) = Perm;
    deriv(1,1) = Perm;
    deriv(2,0) = 1./Perm;
    deriv(3,1) = 1./Perm;
}

template<class TVar>
void TLaplaceExampleSmooth::graduxy(const TPZVec<TVar> &x, TPZVec<TVar> &grad)
{
    TPZManVector<Fad<TVar>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        Fad<TVar> temp = Fad<TVar>(2,i,x[i]);
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<TVar>,3> result(1);
    uxy(xfad,result);
    grad.resize(2);
    for (int i=0; i<2; i++)
    {
        grad[i] = result[0].d(i);
    }
    //    for(int i=0; i<2; i++)
    //    {
    //        std::cout << "result " << result[i] << " dx " << result[i].dx() << std::endl;
    //        for (int j=0; j<2; j++) {
    //            std::cout << "i = " << i << " j = " << j << " grad " << grad(i,j) << " dx " << grad(i,j).dx() << std::endl;
    //        }
    //    }
}

void TLaplaceExampleSmooth::GradU(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu)
{
    TPZManVector<Fad<REAL>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        Fad<REAL> temp = Fad<REAL>(2,i,x[i]);
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<REAL>,3> result(2);
    uxy(xfad,result);
    gradu.Redim(2,1);
    u[0] = result[0].val();
    for (int i=0; i<2; i++) {
        for (int j=0; j<1; j++)
        {
            gradu(i,j) = result[j].d(i);
        }
    }
    
}

template<class TVar>
void TLaplaceExampleSmooth::Sigma(const TPZVec<TVar> &x, TPZVec<TVar> &sigma)
{
    TPZManVector<TVar,3> grad;
    TVar Perm;
    Permeability(x, Perm);
    graduxy(x,grad);
    sigma.resize(2);
    sigma[0] = -Perm*grad[0];
    sigma[1] = -Perm*grad[1];
    
}

template<class TVar>
void TLaplaceExampleSmooth::DivSigma(const TPZVec<TVar> &x, TVar &divsigma)
{
    TPZManVector<Fad<TVar>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        xfad[i] = Fad<TVar>(2,i,x[i]);
    }
    TPZManVector<Fad<TVar>, 3> sigma(2);
    Sigma(xfad,sigma);
    //    for(int i=0; i<2; i++)
    //    {
    //        for (int j=0; j<2; j++) {
    //            std::cout << "i = " << i << " j = " << j <<  " sigma " << sigma(i,j) <<  " " << sigma(i,j).dx() << std::endl;
    //        }
    //    }
    
    divsigma = sigma[0].dx(0)+sigma[1].dx(1);
    
}


TPZAutoPointer<TPZFunction<STATE> > TLaplaceExampleSmooth::ForcingFunction()
{
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Force,5);
    dummy->SetPolynomialOrder(5);
    TPZAutoPointer<TPZFunction<STATE> > result(dummy);
    return result;
}

TPZAutoPointer<TPZFunction<STATE> > TLaplaceExampleSmooth::ValueFunction()
{
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(GradU,5);
    dummy->SetPolynomialOrder(5);
    TPZAutoPointer<TPZFunction<STATE> > result(dummy);
    return result;
    
}

TPZAutoPointer<TPZFunction<STATE> > TLaplaceExampleSmooth::ConstitutiveLawFunction()
{
    TPZAutoPointer<TPZFunction<STATE> > result;
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(PermeabilityDummy,4);
    dummy->SetPolynomialOrder(4);
    result = TPZAutoPointer<TPZFunction<STATE> >(dummy);
    return result;
    
}

#endif

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.meshgen"));
#endif

TPZGeoMesh *MalhaGeomFredQuadrada(int nelx, int nely, TPZVec<REAL> &x0, TPZVec<REAL> &x1, TPZVec<int64_t> &coarseindices, int ndiv)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    int dimension = 2;
    gmesh->SetDimension(dimension);
    TPZManVector<int,2> nx(2,3);
    nx[0] = nelx;
    nx[1] = nely;
    TPZGenGrid gengrid(nx, x0, x1);
    gengrid.SetRefpatternElements(true);
    gengrid.Read(gmesh, 1);
    gengrid.SetBC(gmesh, 4, -1);
    gengrid.SetBC(gmesh, 5, -2);
    gengrid.SetBC(gmesh, 6, -3);
    gengrid.SetBC(gmesh, 7, -4);
    
    // refine a random element
    if(1)
    {
        TPZManVector<TPZGeoEl *,10> gelsub;
        int64_t nel = gmesh->NElements();
        TPZGeoEl *gel = gmesh->Element(rand()%nel);
        while (gel->Dimension() != gmesh->Dimension()) {
            gel = gmesh->Element(rand()%nel);
        }
        gel->Divide(gelsub);
        
    }
    int64_t nel = gmesh->NElements();
    
    coarseindices.resize(nel);
    int64_t elcount = 0;
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->HasSubElement() ||  gel->Dimension() != dimension) {
            continue;
        }
        coarseindices[elcount] = el;
        elcount++;
    }
    coarseindices.resize(elcount);
    
#ifdef PZDEBUG
    {
        std::ofstream file("GMeshFred.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    }
#endif
    if(0)
    {
        
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (!gel->HasSubElement() &&  gel->Type() == EQuadrilateral) {
                TPZManVector<TPZGeoEl *,12> subs;
                gel->Divide(subs);
            }
        }
        nel = gmesh->NElements();
        
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (!gel->HasSubElement() &&  gel->Type() == EOned) {
                TPZAutoPointer<TPZRefPattern> refpat = TPZRefPatternTools::PerfectMatchRefPattern(gel);
                if (!refpat) {
                    DebugStop();
                }
                gel->SetRefPattern(refpat);
                TPZManVector<TPZGeoEl *,12> subs;
                gel->Divide(subs);
            }
        }
    }
    TPZCheckGeom geom(gmesh);
    geom.UniformRefine(ndiv);
    //    InsertInterfaceElements(gmesh,1,2);
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    
    return gmesh;
}

void SolveProblem(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<TPZAutoPointer<TPZCompMesh> > compmeshes, TAnalyticSolution *analytic, std::string prefix, TRunConfig config)
{
    //calculo solution
    bool shouldrenumber = true;
    TPZAnalysis an(cmesh,shouldrenumber);
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh.operator->());
    strmat.SetNumThreads(0);
    
#else
    TPZSkylineStructMatrix strmat(cmesh.operator->());
    strmat.SetNumThreads(0);
#endif
    
    
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    std::cout << "Assembling\n";
    an.Assemble();
    if(0)
    {
        std::string filename = prefix;
        filename += "_Global.nb";
        std::ofstream global(filename.c_str());
        TPZAutoPointer<TPZStructMatrix> strmat = an.StructMatrix();
        an.Solver().Matrix()->Print("Glob = ",global,EMathematicaInput);
        an.Rhs().Print("Rhs = ",global,EMathematicaInput);
    }
    std::cout << "Solving\n";
    an.Solve();
    std::cout << "Finished\n";
    an.LoadSolution(); // compute internal dofs
                       //    an.Solution().Print("sol = ");
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(compmeshes, cmesh);
    
#ifdef PZDEBUG
    {
        std::ofstream out(prefix+"_MeshWithSol.txt");
        cmesh->Print(out);
    }
#endif
    
    //    TPZBuildMultiphysicsMesh::TransferFromMeshes(cmeshes, an.Mesh());
    //    for (int i=0; i<cmeshes.size(); i++) {
    //        cmeshes[i]->Solution().Print("sol = ");
    //    }
    //    cmeshes[0]->Solution().Print("solq = ");
    //    cmeshes[1]->Solution().Print("solp = ");
    std::string plotfile1,plotfile2;
    {
        std::stringstream sout;
        sout << prefix << "Approx-";
        config.ConfigPrint(sout) << "_dim1.vtk";
        plotfile1 = sout.str();
    }
    {
        std::stringstream sout;
        sout << prefix << "Approx-";
        config.ConfigPrint(sout) << "_dim2.vtk";
        plotfile2 = sout.str();
    }
    std::cout << "plotfiles " << plotfile1.c_str() << " " << plotfile2.c_str() << std::endl;
    TPZStack<std::string> scalnames,vecnames;
    TPZMaterial *mat = cmesh->FindMaterial(1);
    if (!mat) {
        DebugStop();
    }
    if (analytic)
    {
        an.SetExact(analytic->Exact());
    }
    if (mat->NStateVariables() == 2)
    {
        scalnames.Push("SigmaX");
        scalnames.Push("SigmaY");
        scalnames.Push("TauXY");
        vecnames.Push("Displacement");
    }
    else if(mat->NStateVariables() == 1)
    {
        scalnames.Push("Pressure");
        vecnames.Push("Flux");
        vecnames.Push("Derivative");
    }
    //    an.DefineGraphMesh(cmesh->Dimension()-1, scalnames, vecnames, plotfile1);
    an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotfile2);
    int resolution = 0;
    //    an.PostProcess(resolution,cmesh->Dimension()-1);
    an.PostProcess(resolution,cmesh->Dimension());
    if(analytic)
    {
        TPZVec<REAL> errors(3,0.);
        an.PostProcessError(errors);
        std::cout << prefix << " - ";
        config.ConfigPrint(std::cout) << " errors computed " << errors << std::endl;
        std::stringstream filename;
        filename << prefix << "Errors.txt";
        std::ofstream out (filename.str(),std::ios::app);
        config.InlinePrint(out);
        out <<  " Energy " << errors[0] << " L2 " << errors[1] << " H1 " << errors[2] << std::endl;
        if (config.newline) {
            out << std::endl;
        }
    }
}

/// Solve the problem composed of a multiphysics mesh composed of compmeshes - applies to MHM and MHM-H(div)
void SolveParabolic(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<TPZAutoPointer<TPZCompMesh> > compmeshes, std::string prefix, TRunConfig config)
{
    
    //calculo solution
    bool shouldrenumber = true;
    TPZAnalysis an(cmesh,shouldrenumber);
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh.operator->());
    strmat.SetNumThreads(0);
    
#else
    TPZSkylineStructMatrix strmat(cmesh.operator->());
    strmat.SetNumThreads(0);
#endif
    
    
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    std::cout << "Assembling\n";
    an.Assemble();
    if(0)
    {
        std::string filename = prefix;
        filename += "_Global.nb";
        std::ofstream global(filename.c_str());
        TPZAutoPointer<TPZStructMatrix> strmat = an.StructMatrix();
        an.Solver().Matrix()->Print("Glob = ",global,EMathematicaInput);
        an.Rhs().Print("Rhs = ",global,EMathematicaInput);
    }
    std::cout << "Solving\n";
    an.Solve();
    std::cout << "Finished\n";
    an.LoadSolution(); // compute internal dofs
                       //    an.Solution().Print("sol = ");
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(compmeshes, cmesh);
    
#ifdef PZDEBUG
    {
        std::ofstream out(prefix+"_MeshWithSol.txt");
        cmesh->Print(out);
    }
#endif
    
    //    TPZBuildMultiphysicsMesh::TransferFromMeshes(cmeshes, an.Mesh());
    //    for (int i=0; i<cmeshes.size(); i++) {
    //        cmeshes[i]->Solution().Print("sol = ");
    //    }
    //    cmeshes[0]->Solution().Print("solq = ");
    //    cmeshes[1]->Solution().Print("solp = ");
    std::string plotfile1,plotfile2;
    {
        std::stringstream sout;
        sout << prefix << "Approx-";
        config.ConfigPrint(sout) << "_dim1.vtk";
        plotfile1 = sout.str();
    }
    {
        std::stringstream sout;
        sout << prefix << "Approx-";
        config.ConfigPrint(sout) << "_dim2.vtk";
        plotfile2 = sout.str();
    }
    std::cout << "plotfiles " << plotfile1.c_str() << " " << plotfile2.c_str() << std::endl;
    TPZStack<std::string> scalnames,vecnames;
    TPZMaterial *mat = cmesh->FindMaterial(1);
    if (!mat) {
        DebugStop();
    }
    if (mat->NStateVariables() == 2)
    {
        scalnames.Push("SigmaX");
        scalnames.Push("SigmaY");
        scalnames.Push("TauXY");
        vecnames.Push("Displacement");
    }
    else if(mat->NStateVariables() == 1)
    {
        scalnames.Push("Pressure");
        vecnames.Push("Flux");
        vecnames.Push("Derivative");
    }
    an.DefineGraphMesh(cmesh->Dimension()-1, scalnames, vecnames, plotfile1);
    an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotfile2);
    int resolution = 0;
    an.PostProcess(resolution,cmesh->Dimension()-1);
    an.SetStep(0);
    an.PostProcess(resolution,cmesh->Dimension());
    
    for (int is=0; is<config.nTimeSteps; is++) {
        an.AssembleResidual();
        an.Solve();
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(compmeshes, cmesh);
        an.PostProcess(resolution,cmesh->Dimension()-1);
        an.SetStep(is+1);
        an.PostProcess(resolution, cmesh->Dimension());

    }
}


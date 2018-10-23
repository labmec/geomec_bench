/*
 *  Monofasico.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 28/07/2017.
 *  Copyright 2017 __MyCompanyName__. All rights reserved.
 *
 */

#include "HidraulicoMonofasicoElastico.h"
#include "pzcheckgeom.h"
#include "pzstack.h"
#include "TPZParSkylineStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzgmesh.h"
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "TPZBndCondWithMem.h"
#include "pzpoisson3d.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmultiphysicselement.h"
#include "TPZInterfaceEl.h"
#include "TPZElasticCriterion.h"
#include "pzporoelastoplasticmem.h"
#include "pzcompelwithmem.h"
#include "../TPZLagrangeInterface.h"
#include "pzelastoplasticanalysis.h"
#include "../TPZSegregatedAnalysisDFN.h"
#include "pzl2projection.h"

#include "../TPZDarcy2DMaterialMem.h"
#include "../TPZInterfaceMemory.h"
#include "../TPZMatFractureBB.h"
#include "../TPZMemoryFracDFN.h"

#include "pzelasmat.h"
#include "pzinterpolationspace.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmultiphysicselement.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzporoelasticmf2d.h"

#include "TPZDarcy2DMaterial.h"
#include "TPZLagrangeMultiplier.h"
#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "pzfmatrix.h"
#include "pzbstrmatrix.h"
#include <TPZGeoElement.h>
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmultiphysicselement.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZGmshReader.h"
#include "../TPZMonoPhasicMemoryDFN.h"
#include "../TPZElastoPlasticMemoryDFN.h"
#include "../TPZMonoPhasicMemoryBCDFN.h"
#include "../TPZElastoPlasticMemoryBCDFN.h"
#include "../TPZPoroElastoPlasticDFN_impl.h"
#include "../TPZMemoryBCDFN.h"
#include "../TPZMemoryDFN.h"

#define TRIANGLEMESH

using namespace std;

const REAL Pi=M_PI;
REAL ftimeatual = 0.;

HidraulicoMonofasicoElastico::HidraulicoMonofasicoElastico()
{
    
    fdim=2; //Dimensão do problema
    fmatID=1; //Materia do elemento volumétrico
    fdimFrac = 1; //Dimensão da fratura
    
    //Materiais das condições de contorno
    fmatBCbott=2;
    fmatBCright=3;
    fmatBCtop=4;
    fmatBCleft=5;
    
    //Número de fraturas do problema:
    fnFrac = 1;
    
    fmatFrac.resize(fnFrac);
    fmatPointLeft.resize(fnFrac);
    fmatPointRight.resize(fnFrac);
    
    for (int i_frac = 0; i_frac < fnFrac; i_frac++) {
        fmatFrac[i_frac] = 6+i_frac;
        fmatPointLeft[i_frac] = 7+3*i_frac;
        fmatPointRight[i_frac] = 8+3*i_frac;
    }
    
    //Material do elemento de interface
    fmatInterface = 500;
    fmatInterfaceLeft = 501;
    fmatInterfaceRight = 502;
    fmatFluxWrap= 503;
    
    //Materiais das condições de contorno (elementos de interface)
    fmatIntBCbott=-11;
    fmatIntBCtop=-12;
    fmatIntBCleft=-13;
    fmatIntBCright=-14;
    
    //Materia de um ponto
    fmatPoint=-5;
    
    //Condições de contorno do problema
    fdirichlet =0;
    fneumann = 1;
    fpenetration=2;
    fpointtype=5;
    fdirichletPress=6;
    fneumdir=10;
    fdirfreey_neum=300;
    fdirneum = 1;
    fmixedneum = 21;
    fmixeddirich = 20;
    fmixedFreeYXdirich = 400;
    fmixedFreeXYdirich = 500;
    
    ftheta=-1;
    fEyoung = 1.;
    fpoisson = 0.;
    falpha = 0.;
    fSe = 0.;
    fperm = 0.;
    fvisc = 0.;
    ffx = 0.;
    ffy = 0.;
    fsign = 0.;
    
    fpref = 0.;
    fLref = 0.;
    fkovervisc = 0.;
    fvalsourceterm = 0.;
    
    finsert_fractures_Q  = true;
    
}

void HidraulicoMonofasicoElastico::Run(int pOrder)
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    //Dados do problema:
    
    HDivPiola = 1;
    
    TPZMaterial::gBigNumber = 1.e12;
    REAL Eyoung = 1.0;
    REAL poisson = 0.0;
    
    REAL rockrho = 0.;
    REAL gravity = 0.;
    REAL fx = 0.;
    REAL fy = gravity*rockrho;
    
    REAL alpha = 0.;
    REAL Se = 0.0;
    REAL perm = 1.;
    REAL visc = 1.;
    REAL sig0 = 0.;
    REAL pini = 0.;
    REAL La = 1.;
    REAL Ly = 5.*La;
    REAL Lx = 8.*La;
    REAL timeT = 0.;
    this->SetParameters(Eyoung, poisson, alpha, Se, perm, visc, fx, fy, sig0);
    finsert_fractures_Q = true;
    
    TPZSimulationData SimData;
    ofstream saidaerro("ErroLoula.txt");

    //Gerando malha geométrica:
    
    TPZGeoMesh *gmesh = CreateGMesh(); //Função para criar a malha geometrica
   
    int n_div = 0;
    
    UniformRef(gmesh,n_div);
    
#ifdef PZDEBUG
    std::ofstream fileg("MalhaGeo.txt"); //Impressão da malha geométrica (formato txt)
    std::ofstream filegvtk("MalhaGeo.vtk"); //Impressão da malha geométrica (formato vtk)
    gmesh->Print(fileg);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk,true);
#endif
    
    //Gerando malha computacional:
    int p = 2;
    TPZCompEl::SetgOrder(p);
    
    TPZSimulationData *simulation_data =  new TPZSimulationData;
    simulation_data->Get_volumetric_material_id().push_back(fmatID);
    if(finsert_fractures_Q){
        for (int i_frac = 0; i_frac < fnFrac; i_frac++) {
            simulation_data->Get_fracture_material_id().push_back(fmatFrac[i_frac]);
        }
    }
    simulation_data->Set_n_threads(0);
    simulation_data->Set_epsilon_res(0.001);
    simulation_data->Set_epsilon_cor(0.001);
    simulation_data->Set_n_iterations(1);
    
    RunningPoroElasticity(gmesh, pOrder, simulation_data);


}

void HidraulicoMonofasicoElastico::RunningPoroElasticity(TPZGeoMesh *gmesh, int pOrder, TPZSimulationData *simulation_data){
    
    //Malha com formulação elástica:
    TPZCompMesh *cmesh_E = CMesh_E(gmesh, pOrder,simulation_data);
    
    if (finsert_fractures_Q) {
        
        std::ofstream filecEbf("MalhaC_E_before.txt"); //Impressão da malha computacional da velocidade (formato txt)
        cmesh_E->Print(filecEbf);
        std::vector<int> fracture_ids;
        for (int i_frac=0; i_frac< fnFrac; i_frac++) {
            fracture_ids.push_back(fmatFrac[i_frac]);
        }
        BreakH1Connectivity(*cmesh_E, fracture_ids); // Insert new connects to represent normal fluxes

        std::ofstream fileg1("MalhaGeo2.txt"); //Impressão da malha geométrica (formato txt)
        gmesh->Print(fileg1);
        std::ofstream filegvtk2("MalhaGeo2.vtk"); //Impressão da malha geométrica (formato vtk)
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk2,true);
        
    }
    {
        std::ofstream filecE("MalhaC_E.txt"); //Impressão da malha computacional da velocidade (formato txt)
        cmesh_E->Print(filecE);
    }

    //Malha com formulação mista (Darcy):
    TPZManVector<TPZCompMesh* , 2 > mesh_vector(2);
    TPZCompMesh *cmesh_M = CMesh_M(mesh_vector, gmesh, pOrder-1, simulation_data);
    
    // ******* write a procedure to set equalize the integration order of the meshes
    cmesh_M->LoadReferences();
    long nel = cmesh_M->NElements();
    
    TPZVec<long> indices;
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh_M->Element(el);
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfcel)
        {
            continue;
        }
        mfcel->InitializeIntegrationRule();
        mfcel->PrepareIntPtIndices();
    }
//    cmesh_M->InitializeBlock();
    
    //Add multiphysics interfaces
    if (finsert_fractures_Q) {
        for (int i_frac = 0; i_frac < fnFrac; i_frac++) {
            AddMultiphysicsInterfaces(*cmesh_M,fmatFrac[i_frac]);
        }
    }
    
    
    {
        std::ofstream filecMbf("MalhaC_M.txt");
        cmesh_M->Print(filecMbf);
    }
    //Cria Análise das malhas
    
    std::string plotfileM("Benchmark_HME_Darcy.vtk");
    std::string plotfileE("Benchmark_HME_Elastoplast.vtk");
    
    TPZStack<std::string> var_names_darcy;
    var_names_darcy.Push("P");
    var_names_darcy.Push("Vx");
    var_names_darcy.Push("Vy");
    
    TPZStack<std::string> var_names_elastoplast;
    var_names_elastoplast.Push("sxx");
    var_names_elastoplast.Push("syy");
    var_names_elastoplast.Push("ux");
    var_names_elastoplast.Push("uy");
    
    TPZSegregatedAnalysisDFN * segregated_analysis = new TPZSegregatedAnalysisDFN;
    
    segregated_analysis->ConfigurateAnalysis(ELDLt, ELU, simulation_data , cmesh_E, cmesh_M, mesh_vector, var_names_elastoplast, var_names_darcy);
    //Depois adicionar variaveis vetoriais Darcy
    
    segregated_analysis->ExecuteTimeEvolution();

    
    if(finsert_fractures_Q) {
        std::cout << "Postprocessing fracture" << std::endl;
        Plot_over_fractures(cmesh_M);
    }
    
    
}


void HidraulicoMonofasicoElastico::UniformRef(TPZGeoMesh * gmesh, int n_div){
    for ( int ref = 0; ref < n_div; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        int n = gmesh->NElements();
        for ( int i = 0; i < n; i++ ){
            TPZGeoEl * gel = gmesh->Element(i);
            gel->Divide (filhos);
        }//for i
    }//ref
}


HidraulicoMonofasicoElastico::~HidraulicoMonofasicoElastico()
{
    
}

void HidraulicoMonofasicoElastico::F_source(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE>& gradu){
    
    f.resize(1);
    const REAL Pi=M_PI;
    
    REAL xv = x[0];
    REAL yv = x[1];
    
    STATE f_x = 0.;
    
    f[0] = f_x;
    
}

TPZAutoPointer <TPZMatrix<STATE> > HidraulicoMonofasicoElastico::MassMatrix(TPZPoroElasticMF2d * mymaterial, TPZCompMesh* mphysics, int nthreads){
    
    mymaterial->SetLastState();
    //TPZSkylineStructMatrix matsp(mphysics);
    //TPZSkylineNSymStructMatrix matsp(mphysics);
    TPZSpStructMatrix matsp(mphysics);
    matsp.SetNumThreads(nthreads);
    
    std::set< int > materialid;
    int matid = mymaterial->MatId();
    materialid.insert(matid);
    matsp.SetMaterialIds (materialid);
    TPZAutoPointer<TPZGuiInterface> guiInterface;
    TPZFMatrix<STATE> Un;
    
    
    //TPZMatrix<REAL> *matK2 = matsp.CreateAssemble(Un,guiInterface);
    
    TPZAutoPointer <TPZMatrix<STATE> > matK2 = matsp.CreateAssemble(Un,guiInterface);
    
    return matK2;
}

void HidraulicoMonofasicoElastico::StiffMatrixLoadVec(TPZPoroElasticMF2d *mymaterial, TPZCompMesh* mphysics, TPZAnalysis &an, TPZFMatrix<STATE> &matK1, TPZFMatrix<STATE> &fvec, int nthreads) {
    
    mymaterial->SetCurrentState();
    //TPZFStructMatrix matsk(mphysics);
    TPZSkylineStructMatrix matsk(mphysics);
    matsk.SetNumThreads(nthreads);
    an.SetStructuralMatrix(matsk);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    //step.SetDirect(ELU);
    an.SetSolver(step);
    an.Run();
    
    matK1 = an.StructMatrix();
    fvec = an.Rhs();
    
}

void HidraulicoMonofasicoElastico::SetParameters(REAL mod_young, REAL mod_poisson, REAL coef_alpha, REAL coef_Se, REAL permeabil_fluido, REAL visc_fluido, REAL fx, REAL fy,REAL sign){
    
    fEyoung = mod_young;
    fpoisson= mod_poisson;
    falpha = coef_alpha;
    fSe = coef_Se;
    fperm = permeabil_fluido;
    fvisc = visc_fluido;
    ffx = fx;
    ffy = fy;
    fsign = sign;
}

void HidraulicoMonofasicoElastico::PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile)
{
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    
    TPZManVector<std::string,10> scalnames(3), vecnames(2);
    //scalnames[0] = "DisplacementX";
    //scalnames[1] = "DisplacementY";
    vecnames[0] = "Displacement";
    scalnames[0] = "SigmaX";
    scalnames[1] = "SigmaY";
    //scalnames[0] = "PorePressure";
    //scalnames[2] = "FluxoY";
    //vecnames[0] = "Fluxo";
    //vecnames[1] = "MinusKMuGradP";
    
    scalnames[2] = "ExactPressure";
    //scalnames[6] = "FluxoX";
    
    //scalnames[4] = "ExactDisplacementX";
    //scalnames[5] = "ExactDisplacementY";
    // scalnames[8] = "ExactSigmaX";
    //scalnames[6] = "ExactSigmaY";
    vecnames[1]  = "ExactFluxo";
    //vecnames[1]  = "ExactDisplacement";
    //vecnames[4] = "MinusKMuGradP";
    
    const int dim = 2;
    int div = 0;
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    //    std::ofstream out("malha.txt");
    //    an.Print("nothing",out);
}

TPZGeoMesh *HidraulicoMonofasicoElastico::CreateGMesh()
{
    int64_t id, index;
    
    //Criando malha geométrica, nós e elementos.
    //Inserindo nós e elementos no objeto malha:
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    //Aqui é implementado um método para malhas criadas no GMSH
    
    //std::string dirname = PZSOURCEDIR;
    std::string grid;
    
    grid = "/Users/pablocarvalho/Documents/GitHub/geomec_bench/Fase_1/Benchmark0a/gmsh/GeometryBench.msh";

    TPZGmshReader Geometry;
    REAL s = 1.0;
    Geometry.SetfDimensionlessL(s);
    Geometry.fPZMaterialId[0]["PointLeft"] = 0;
    Geometry.fPZMaterialId[0]["PointRight"] = 0;
    Geometry.fPZMaterialId[1]["bottom"] = fmatBCbott;
    Geometry.fPZMaterialId[1]["right"] = fmatBCright;
    Geometry.fPZMaterialId[1]["top"] = fmatBCtop;
    Geometry.fPZMaterialId[1]["left"] = fmatBCleft;
    if (finsert_fractures_Q) {
        Geometry.fPZMaterialId[1]["frac"] = fmatFrac[0];
    }
    Geometry.fPZMaterialId[2]["Omega"] = fmatID;
    gmesh = Geometry.GeometricGmshMesh(grid);
    
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();
    
    gmesh->BuildConnectivity();
    
    //    int n_div = 0;
    //    UniformRefine(gmesh,n_div);
    ofstream bf("before.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
    return gmesh;

}

void HidraulicoMonofasicoElastico::Sol_exact(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &deriv){

    REAL y = ptx[1];
    sol.Resize(2, 0.);// ux, uy;
    
    sol[1] = 0.;
    deriv.Resize(2,2);//sigx, sigxy, sigyx, sigy
    deriv(0,0) = deriv(0,1) = deriv(1,0) = deriv(1,1) = 0.;
    
    
}

TPZCompMesh *HidraulicoMonofasicoElastico::CMesh_E(TPZGeoMesh *gmesh, int pOrder, TPZSimulationData *sim_data)
{
    
    // Criando malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    if(pOrder<2) DebugStop();
    
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(fdim);
    cmesh->SetAllCreateFunctionsContinuousWithMem();
    
    // 1 - Material volumétrico
    int nstate = 2;
    TPZVec<REAL> force(fdim,0.);
    //REAL E = 0;
    //REAL poisson = 0;
    
    TPZPoroElastoPlasticDFN<TPZElasticCriterion,TPZMemoryDFN> *material;
    //material = new TPZElasticityMaterial(fmatID, fEyoung, fpoisson, ffx, ffy, planestress);
    material = new TPZPoroElastoPlasticDFN<TPZElasticCriterion , TPZMemoryDFN> (fmatID,fdim);
    TPZElasticCriterion obj;
    
    TPZElasticResponse er;
    er.SetUp(fEyoung,fpoisson);
    obj.SetElasticResponse(er);
    material->SetPlasticity(obj);
    material->SetAlpha(1.);
    material->SetSimulationData(sim_data);
    cmesh->InsertMaterialObject(material);
    
    int null_dir_dirichlet = 3;
    
    // 1 - Condições de contorno
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    val2(1,0)= 1;
    
    TPZBndCondWithMem<TPZMemoryBCDFN> * BCond1 = new TPZBndCondWithMem<TPZMemoryBCDFN>(material, fmatBCbott, null_dir_dirichlet, val1, val2);
    val2.Zero();
    val2(1,0)= -1.;
    TPZBndCondWithMem<TPZMemoryBCDFN> * BCond2 = new TPZBndCondWithMem<TPZMemoryBCDFN>(material, fmatBCtop, fneumann, val1, val2);
    val2.Zero();
    val2(0,0)=1;
    TPZBndCondWithMem<TPZMemoryBCDFN> * BCond3 = new TPZBndCondWithMem<TPZMemoryBCDFN>(material, fmatBCright, null_dir_dirichlet, val1, val2);
    val2.Zero();
    val2(0,0)=1;
    TPZBndCondWithMem<TPZMemoryBCDFN> * BCond4 = new TPZBndCondWithMem<TPZMemoryBCDFN>(material, fmatBCleft, null_dir_dirichlet, val1, val2);
    
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    
    // 2 - Material Fraturas
    if (finsert_fractures_Q) {
     
        TPZMatFractureBB<TPZMemoryFracDFN> *materialFrac;
        for (int i_frac = 0; i_frac < fnFrac; i_frac++) {
            materialFrac = new TPZMatFractureBB<TPZMemoryFracDFN>(fmatFrac[i_frac],fdimFrac,nstate);
            cmesh->InsertMaterialObject(materialFrac);
        }
        
        // 2 - Material Lagrange nas interfaces
        TPZLagrangeInterface<TPZInterfaceMemory> *matInterLeft = new TPZLagrangeInterface<TPZInterfaceMemory>(fmatInterfaceLeft, fdimFrac, nstate);
        matInterLeft->SetMultiplier(-1);
        matInterLeft->SetSimulationData(sim_data);
        cmesh->InsertMaterialObject(matInterLeft);
        
        TPZLagrangeInterface<TPZInterfaceMemory> *matInterRight = new TPZLagrangeInterface<TPZInterfaceMemory>(fmatInterfaceRight, fdimFrac, nstate);
        matInterRight->SetMultiplier(1);
        matInterRight->SetSimulationData(sim_data);
        cmesh->InsertMaterialObject(matInterRight);
        
        // 2 - Condições de contorno
        //    for (int i_frac = 0; i_frac < fnFrac; i_frac++) {
        //        //Mat Frac:
        //        TPZMaterial * bc_frac_right = materialFrac->CreateBC(materialFrac, fmatPointRight[i_frac] , fdirichlet, val1, val2);
        //        cmesh->InsertMaterialObject(bc_frac_right);
        //        TPZMaterial * bc_frac_left = materialFrac->CreateBC(materialFrac, fmatPointLeft[i_frac] , fdirichlet, val1, val2);
        //        cmesh->InsertMaterialObject(bc_frac_left);
        //    }
        
        // 2 - Criando material para FluxWrap
        //    TPZBndCond * bc_fracture_wrap;
        //    bc_fracture_wrap = material->CreateBC(material,fmatFluxWrap,fdirichlet,val1,val2);
        //    cmesh->InsertMaterialObject(bc_fracture_wrap);
        
    }
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    
    std::set<int> matids;
    matids.insert(fmatID);
    matids.insert(fmatBCbott);
    matids.insert(fmatBCright);
    matids.insert(fmatBCtop);
    matids.insert(fmatBCleft);


    if (finsert_fractures_Q) {
        for (int i_frac = 0; i_frac < fnFrac; i_frac++) {
            matids.insert(fmatFrac[i_frac]);
        }
        matids.insert(fmatInterfaceRight);
        matids.insert(fmatInterfaceLeft);
    }
    
    cmesh->ApproxSpace().CreateWithMemory(true);// Force the creation of interfaces with memory.
    cmesh->AutoBuild(matids);
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->InitializeBlock();
    
    return cmesh;
}

TPZCompMesh *HidraulicoMonofasicoElastico::CMesh_q(TPZGeoMesh *gmesh, int pOrder)
{
    
    // Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    int nstate = 1;
    TPZVec<STATE> sol;
    cmesh->SetDefaultOrder(pOrder);//Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim);//Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsHDiv(); //Criando funções HDIV:
    TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xbin(1,1,0.), xfin(1,1,0.);
    
    // 1 - Material volumétrico
    TPZL2Projection *material = new TPZL2Projection(fmatID,fdim,nstate,sol);
    cmesh->InsertMaterialObject(material);
    
    // 1 - Condições de contorno
    TPZFMatrix<STATE> val1(1,1,0.), val2(3,1,0.);
    TPZMaterial * BCond0 = material->CreateBC(material, fmatBCbott, fdirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCond0);
    
    TPZMaterial * BCond1 = material->CreateBC(material, fmatBCtop, fdirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCond1);
    
    TPZMaterial * BCond2 = material->CreateBC(material, fmatBCright, fdirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCond2);
    
    TPZMaterial * BCond3 = material->CreateBC(material, fmatBCleft, fdirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCond3);
    
    // 2 - Material Fraturas
    if (finsert_fractures_Q) {
        TPZMat1dLin *materialFrac;
        for (int i_frac = 0; i_frac < fnFrac; i_frac++) {
            materialFrac = new TPZMat1dLin(fmatFrac[i_frac]);

            materialFrac->SetMaterial(xkin, xcin, xbin, xfin);
            cmesh->InsertMaterialObject(materialFrac);

        // 2 - Condições de contorno
         TPZMaterial * BCond4 = materialFrac->CreateBC(materialFrac, fmatPointRight[i_frac] , fdirichlet, val1, val2);
         cmesh->InsertMaterialObject(BCond4);
         TPZMaterial * BCond5 = materialFrac->CreateBC(materialFrac, fmatPointLeft[i_frac] , fdirichlet, val1, val2);
         cmesh->InsertMaterialObject(BCond5);

        }
        
        // 2 - Criando material para FluxWrap
        TPZBndCond *FluxWrapBC;
        FluxWrapBC = material->CreateBC(material,fmatFluxWrap,fdirichlet,val1,val2);
        cmesh->InsertMaterialObject(FluxWrapBC);
        
    }
    

    // Criando elementos computacionais que gerenciarão o espaco de aproximacao da malha:
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
    }
    
    std::set<int> matids;
    matids.insert(fmatID);
    matids.insert(fmatBCbott);
    matids.insert(fmatBCright);
    matids.insert(fmatBCtop);
    matids.insert(fmatBCleft);
    
    if (finsert_fractures_Q) {
        matids.insert(fmatFluxWrap);
    }
    
    cmesh->AutoBuild(matids);

    if (finsert_fractures_Q) {
        gmesh->ResetReference();
        matids.clear();
        
        for (int i_frac = 0; i_frac < fnFrac; i_frac++) {
            matids.insert(fmatFrac[i_frac]);
            matids.insert(fmatPointLeft[i_frac]);
            matids.insert(fmatPointRight[i_frac]);
        }
    
        cmesh->SetDimModel(fdimFrac);
        cmesh->SetAllCreateFunctionsHDiv();
        cmesh->AutoBuild(matids);
    }
    
    cmesh->AdjustBoundaryElements();
    cmesh->ExpandSolution();
    cmesh->SetDimModel(fdim);
    cmesh->SetAllCreateFunctionsHDiv();

    return cmesh;
    
}


TPZCompMesh *HidraulicoMonofasicoElastico::CMesh_p(TPZGeoMesh *gmesh, int pOrder)
{
    
    // @omar::
    
    //pOrder--; // Space restriction apapapa
    
    // Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    int nstate = 1;
    TPZVec<STATE> sol;
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
    
    // 1 - Material volumétrico
    TPZL2Projection *material = new TPZL2Projection(fmatID,fdim,nstate,sol);//criando material que implementa a formulacao fraca do problema modelo
    cmesh->InsertMaterialObject(material); //Insere material na malha
    
    // 2 - Material Fraturas
    if (finsert_fractures_Q) {
        TPZMat2dLin *materialFrac;
        for (int i_frac = 0; i_frac < fnFrac; i_frac++) {
            materialFrac = new TPZMat2dLin(fmatFrac[i_frac]);
            cmesh->InsertMaterialObject(materialFrac);
        }
        materialFrac->SetMaterial(xkin, xcin, xfin);
    }
    
    // Criando elementos computacionais que gerenciarão o espaco de aproximacao da malha:

    std::set<int> matids;
    matids.insert(fmatID);
    
    cmesh->SetDimModel(fdim);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    if (finsert_fractures_Q) {
        //gmesh->ResetReference();
        //matids.clear();
        for (int i_frac = 0; i_frac < fnFrac; i_frac++) {
            matids.insert(fmatFrac[i_frac]);
        }
        cmesh->SetDimModel(fdimFrac); //Insere dimensão do modelo
    }
    
    cmesh->AutoBuild(matids);
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    return cmesh;

}

TPZCompMesh *HidraulicoMonofasicoElastico::CMesh_M(TPZManVector<TPZCompMesh* , 2 >  &mesh_vector, TPZGeoMesh *gmesh, int pOrder, TPZSimulationData *sim_data){
    
    int nstate = 2;
    
    TPZCompMesh *cmesh_q = CMesh_q(gmesh,pOrder);
    TPZCompMesh *cmesh_p = CMesh_p(gmesh,pOrder);
    
    {
        std::ofstream filecq("MalhaC_q.txt"); //Impressão da malha computacional da velocidade (formato txt)
        std::ofstream filecp("MalhaC_p.txt"); //Impressão da malha computacional da pressão (formato txt)
        cmesh_q->Print(filecq);
        cmesh_p->Print(filecp);
    }

    if (finsert_fractures_Q) {
        for (int i_frac=0; i_frac< fnFrac; i_frac++) {
            BreakConnectivity(*cmesh_q, fmatFrac[i_frac]); // Insert new connects to represent normal fluxes
        }
    }
    
    // Criando malha computacional:
    int bc_inte_order = 10;
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim); //Insere dimensão do modelo
    
    // 1 - Material volumétrico
    TPZDarcy2DMaterialMem<TPZMemoryDFN> *material = new TPZDarcy2DMaterialMem<TPZMemoryDFN> (fmatID,fdim,1,1);
    TPZFMatrix<REAL> K(fdim,fdim),invK(fdim,fdim);
    K.Zero();
    invK.Zero();
    
    K(0,0)=3.38801e-7;
    K(1,1)=2.566e-11;
    invK(0,0)=1./K(0,0);
    invK(1,1)=1./K(1,1);

    TPZAutoPointer<TPZFunction<STATE> > solp = new TPZDummyFunction<STATE> (Sol_exact,5);
    material->SetForcingFunctionExact(solp);
    material->SetSimulationData(sim_data);
    cmesh->InsertMaterialObject(material);

    // 1 - Condições de contorno
    TPZFMatrix<STATE> val1(1,1,0.), val2(3,1,0.);
    STATE Pjusante = 54.9;
    STATE Pmontante = 55.0;
    
    val1(0,0) = 0.; //botton
    TPZBndCondWithMem<TPZMemoryBCDFN> * BCond0 = new TPZBndCondWithMem<TPZMemoryBCDFN>(material, fmatBCbott, fneumann, val1, val2);
    cmesh->InsertMaterialObject(BCond0);
    
    val1(0,0) = 0.; //top
    TPZBndCondWithMem<TPZMemoryBCDFN> * BCond1 = new TPZBndCondWithMem<TPZMemoryBCDFN>(material, fmatBCtop, fneumann, val1, val2);
    cmesh->InsertMaterialObject(BCond1);
    
    val1(0,0) = Pjusante; // right
    TPZBndCondWithMem<TPZMemoryBCDFN> * BCond2 = new TPZBndCondWithMem<TPZMemoryBCDFN>(material, fmatBCright, fdirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCond2);
    
    val1(0,0) = Pmontante; // left
    TPZBndCondWithMem<TPZMemoryBCDFN> * BCond3 = new TPZBndCondWithMem<TPZMemoryBCDFN>(material, fmatBCleft, fdirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCond3);
    val1(0,0) = 0.0;
    
    if (finsert_fractures_Q) {
        // 2 - Material Fraturas
        TPZDarcy2DMaterialMem<TPZMemoryFracDFN> *materialFrac;
        for (int i_frac = 0; i_frac < fnFrac; i_frac++) {
            materialFrac = new TPZDarcy2DMaterialMem<TPZMemoryFracDFN> (fmatFrac[i_frac],fdimFrac,1,1);
            REAL kf = 4.68789e-4;
            REAL Dyf = 6.5e-5;
            materialFrac->SetPermeability(kf*Dyf);
            materialFrac->SetSimulationData(sim_data);
            cmesh->InsertMaterialObject(materialFrac);
            
            // 2 - Condições de contorno
            val1(0,0) =  Pjusante; // right
            TPZBndCondWithMem<TPZMemoryBCDFN> * BCond4 = new TPZBndCondWithMem<TPZMemoryBCDFN>(materialFrac, fmatPointRight[i_frac], fdirichlet, val1, val2);
            cmesh->InsertMaterialObject(BCond4);

            val1(0,0) = Pmontante; // left
            TPZBndCondWithMem<TPZMemoryBCDFN> * BCond5 = new TPZBndCondWithMem<TPZMemoryBCDFN>(materialFrac, fmatPointLeft[i_frac], fdirichlet, val1, val2);
            cmesh->InsertMaterialObject(BCond5);
        }
        
        // 2 - Material Lagrange nas interfaces
        
        TPZLagrangeInterface<TPZInterfaceMemory> *MatLagrange = new TPZLagrangeInterface<TPZInterfaceMemory>(fmatInterface,fdimFrac,1);
        MatLagrange->SetSimulationData(sim_data);
        cmesh->InsertMaterialObject(MatLagrange);
        
        // 2 - Flux Warap
        val1.Zero();
        val1.Zero();
        TPZBndCond *FluxWrapBC = material->CreateBC(material,fmatFluxWrap,fdirichlet,val1,val2);
        cmesh->InsertMaterialObject(FluxWrapBC);
    }
    

    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha:
    cmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    cmesh->ApproxSpace().CreateWithMemory(true);

    std::set<int> matids;
    matids.insert(fmatID);
    matids.insert(fmatBCbott);
    matids.insert(fmatBCright);
    matids.insert(fmatBCtop);
    matids.insert(fmatBCleft);
    
    if (finsert_fractures_Q) {
        for (int i_frac = 0; i_frac < fnFrac; i_frac++) {
            matids.insert(fmatFrac[i_frac]);
            matids.insert(fmatPointLeft[i_frac]);
            matids.insert(fmatPointRight[i_frac]);
        }
    }
    matids.insert(fmatInterface);
    
    cmesh->AutoBuild(matids);

    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->ApproxSpace().CreateWithMemory(false);
    cmesh->LoadReferences();
    cmesh->AutoBuild();
    
    cmesh->CleanUpUnconnectedNodes();
    
    mesh_vector[0] = cmesh_q;
    mesh_vector[1] = cmesh_p;
    
    TPZBuildMultiphysicsMesh::AddElements(mesh_vector, cmesh);
    TPZBuildMultiphysicsMesh::AddConnects(mesh_vector, cmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(mesh_vector, cmesh);
    
    return cmesh;
    
}


void HidraulicoMonofasicoElastico::Plot_over_fractures(TPZCompMesh *cmesh){
    // Counting for n_fractures
    int ncel = cmesh->NElements();
    TPZStack<int> frac_indexes;
    for (int icel = 0; icel<ncel; icel++) {
        TPZCompEl *cel=cmesh->Element(icel);
        if (!cel) {
            DebugStop();
        }
        TPZGeoEl *gel=cel->Reference();
        if (!gel) {
            DebugStop();
        }
        if (gel->MaterialId()!=fmatFrac[0]) {
            continue;
        }
        frac_indexes.Push(cel->Index());
    }
    
    int n_frac_cels = frac_indexes.size();
    int var_p = 0;
    TPZVec<REAL> par_xi(1,0.0);
    TPZManVector<STATE,3> sol,x(3,0.0);
    TPZFMatrix<REAL> pressure(n_frac_cels*3,3,0.0);
    
    TPZManVector<REAL,3> par_vals(3,0.0);
    par_vals[0] = -1.0;
    par_vals[1] =  0.0;
    par_vals[2] = +1.0;
    
    
    for (int ifrac = 0; ifrac < n_frac_cels; ifrac++) {
        
        TPZCompEl *cel= cmesh->Element(frac_indexes[ifrac]);
        TPZGeoEl *gel=cel->Reference();
        for (int ip = 0; ip < par_vals.size(); ip++) {
            par_xi[0] = par_vals[ip];
            gel->X(par_xi, x);
            cel->Solution(par_xi, var_p, sol);
            pressure(ifrac*3+ip,0) = x[0];
            pressure(ifrac*3+ip,1) = x[1];
            pressure(ifrac*3+ip,2) = sol[0];
        }
        
    }
    
    pressure.Print("pf = ",std::cout,EMathematicaInput);
}


void HidraulicoMonofasicoElastico::BreakConnectivity(TPZCompMesh &cmesh, int matId)
{
    std::set<int> boundaries_ids;
    boundaries_ids.insert(fmatBCbott);
    boundaries_ids.insert(fmatBCleft);
    boundaries_ids.insert(fmatBCtop);
    boundaries_ids.insert(fmatBCright);
    
    TPZFractureInsertion fracture(cmesh.Reference(),matId,boundaries_ids);
    fracture.OpenFractureOnHdiv(&cmesh,fmatFluxWrap);
}

void HidraulicoMonofasicoElastico::BreakH1Connectivity(TPZCompMesh &cmesh, std::vector<int> fracture_ids)
{
    std::set<int> boundaries_ids;
    boundaries_ids.insert(fmatBCbott);
    boundaries_ids.insert(fmatBCleft);
    boundaries_ids.insert(fmatBCtop);
    boundaries_ids.insert(fmatBCright);
    
    for (unsigned int i_f = 0; i_f <  fracture_ids.size(); i_f++) {
        TPZFractureInsertion fracture(cmesh.Reference(),fracture_ids[i_f],boundaries_ids);
        fracture.ClassifyNeighboursofPivots();
        fracture.OpenFractureOnH1(&cmesh); // (ok)
        fracture.SetDiscontinuosFrac(&cmesh); // (ok)
        fracture.SetInterfaces(&cmesh, fmatInterfaceLeft, fmatInterfaceRight);
    }
    cmesh.ComputeNodElCon();
    
}

void HidraulicoMonofasicoElastico::AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int mat_Frac)
{

    TPZGeoMesh *gmesh = cmesh.Reference();
    std::set<int> velmatid;
    velmatid.insert(mat_Frac);
    
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        int matid = gel->MaterialId();
        if(velmatid.find(matid) == velmatid.end())
        {
            continue;
        }
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        TPZGeoElSide neighbour = gelside.Neighbour();
        while (neighbour != gelside) {
            
            TPZManVector<int64_t,3> LeftElIndices(1,0.),RightElIndices(1,0.);
            LeftElIndices[0]=0;
            RightElIndices[0]=1;
            
            if (neighbour.Element()->Dimension() == 1 && neighbour.Element()->MaterialId() == fmatFluxWrap) { //oioioi IDFlux -> ID
                // create an interface element
                TPZCompElSide celside = gelside.Reference();
                TPZCompElSide Wrapneigh = neighbour.Reference();
                if (!celside || !Wrapneigh) {
                    DebugStop();
                }
                std::cout << "Created an element between volumetric element " << neighbour.Element()->Index() <<
                " side " << neighbour.Side() <<
                " and interface element " << gelside.Element()->Index() << std::endl;
                TPZGeoElBC gelbc(gelside,fmatInterface);
                int64_t index;
                TPZMultiphysicsInterfaceElement *intf = new
                TPZMultiphysicsInterfaceElement(cmesh,gelbc.CreatedElement(),index,Wrapneigh,celside);
                intf->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
                
            }
            neighbour = neighbour.Neighbour();
        }
    }
    
}

void HidraulicoMonofasicoElastico::AdjustIntegrationOrder(TPZSimulationData * sim_data, TPZCompMesh * cmesh_elastoplast, TPZCompMesh * cmesh_darcy){
    
    int nel_elastoplast = cmesh_elastoplast->NElements();
    int nel_darcy = cmesh_darcy->NElements();
    
    if (nel_elastoplast != nel_darcy) {
        std::cout << "The geometrical partitions are not the same." << std::endl;
        DebugStop();
    }
    
    int p_order_elastoplast = sim_data->Get_elasticity_order();
    int p_order_darcy = sim_data->Get_darcy_order();
    
    int int_order = 0;
    
//    if (sim_data->Get_is_dual_formulation_Q()) {
//        if (p_order_darcy > p_order_elastoplast) {
//            int_order = 2*(p_order_darcy);
//        }else{
//            int_order = 2*(p_order_elastoplast);
//        }
//
//        if(p_order_darcy == p_order_elastoplast){
//            int_order = 2*(p_order_darcy+1);
//        }
//    }else{
        if (p_order_darcy > p_order_elastoplast) {
            int_order = 2*(p_order_darcy+1);
        }else{
            int_order = 2*(p_order_elastoplast);
        }
        
        if(p_order_darcy == p_order_elastoplast){
            int_order = 2*(p_order_darcy);
        }
//    }
    
    std::vector<TPZIntPoints *> int_rule_vec;
    for (long el = 0; el<nel_elastoplast; el++) {
        TPZCompEl *cel = cmesh_elastoplast->Element(el);
        if (!cel) {
            continue;
        }
        cel->SetIntegrationRule(int_order);
        cel->PrepareIntPtIndices();
        const TPZIntPoints & rule = cel->GetIntegrationRule();
        TPZIntPoints * rule_copy = rule.Clone();
        int_rule_vec.push_back(rule_copy);
    }
    
#ifdef PZDEBUG
    std::ofstream out_elastoplast("Cmesh_Geomechanics_adjusted.txt");
    cmesh_elastoplast->Print(out_elastoplast);
#endif
    
    int counter = 0;
    for (long el = 0; el < nel_darcy; el++) {
        TPZCompEl *cel = cmesh_darcy->Element(el);
//        if (sim_data->Get_is_dual_formulation_Q()) {
//            TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
//            if (!mfcel) {
//                continue;
//            }
//            mfcel->SetIntegrationRule(int_rule_vec[counter]);
//            mfcel->PrepareIntPtIndices();
//        }else{
            cel->SetIntegrationRule(int_order);
//        }
        counter++;
    }
    
//    if (sim_data->Get_is_dual_formulation_Q()) {
//        cmesh_reservoir->CleanUpUnconnectedNodes();
//    }
    
#ifdef PZDEBUG
    std::ofstream out_darcy("CmeshReservoir_adjusted.txt");
    cmesh_darcy->Print(out_darcy);
#endif
    
    
}




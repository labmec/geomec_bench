/*
 *  HidraulicoMonofasicoElastico.h
 *  PZ
 *
 *  Created by Pablo Carvalho on 28/07/2017.
 *  Copyright 2017 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __PZ__HidraulicoMonofasicoElastico__
#define __PZ__HidraulicoMonofasicoElastico__

#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "HidraulicoMonofasicoElastico.h"
#include "pzporoelasticmf2d.h"

#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "pzfmatrix.h"
#include "pzbstrmatrix.h"
#include <TPZGeoElement.h>
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmat2dlin.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzanalysis.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzpostprocanalysis.h"
#include "../TPZFractureInsertion.h"
#include "../TPZSimulationData.h"

using namespace std;
using namespace pzshape;

class HidraulicoMonofasicoElastico{
private:
    
    int fdim; //Dimensão do problema
    int fmatID; //Materia do elemento volumétrico
    int fdimFrac; //Dimensão da fratura
    
    //Materiais das condições de contorno
    int fmatBCbott;
    int fmatBCtop;
    int fmatBCleft;
    int fmatBCright;
    vector<int> fmatFrac;
    vector<int> fmatPointLeft;
    vector<int> fmatPointRight;
    
    std::map<REAL, REAL> fFracOrient;
    int fnFrac;
    
    //int fmatFrac;
    //int fmatPointLeft;
    //int fmatPointRight;
    
    //Material do elemento de interface

    vector<int>  fmatInterface;
    vector<int>  fmatInterfaceLeft;
    vector<int>  fmatInterfaceRight;
    vector<int>  fmatFluxWrap;
    
    //Materiais das condições de contorno (elementos de interface)
    int fmatIntBCbott;
    int fmatIntBCtop;
    int fmatIntBCleft;
    int fmatIntBCright;
    
    //Materia de um ponto
    int fmatPoint;
    
    //Condições de contorno do problema
    int fdirichlet = 0;
    int fneumann = 1;
    int fmixed = 2;
    int fpenetration;
    int fpointtype;
    int fdirichletPress;
    
    int  fneumdir;
    int  fdirfreey_neum;
    int  fdirneum;
    int  fmixedneum;
    int  fmixeddirich;
    
    int fmixedFreeXYdirich;
    int fmixedFreeYXdirich;
    
    REAL fEyoung;
    REAL fpoisson;
    REAL falpha;
    REAL fSe;
    REAL fperm;
    REAL fvisc;
    REAL ffx;
    REAL ffy;
    REAL fsign;
    
    REAL fpref;
    REAL fLref;
    REAL fkovervisc;
    
    REAL fvalsourceterm;
    
    int ftheta;
    
    bool finsert_fractures_Q = true;
    
    TPZVec<TPZFractureInsertion > fractureInsert;
    
    
public:

    HidraulicoMonofasicoElastico();
    
    void Run(int pOrder);
    
    void RunningPoroElasticity(TPZGeoMesh *gmesh, int pOrder, TPZSimulationData *simulation_data);
    
    ~HidraulicoMonofasicoElastico();
    
    static void FunctionStress(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE>& gradu);
    
    void StiffMatrixLoadVec(TPZPoroElasticMF2d *mymaterial, TPZCompMesh* mphysics, TPZAnalysis &an, TPZFMatrix<STATE> &matK1, TPZFMatrix<STATE> &fvec, int nthreads);
    
    TPZAutoPointer <TPZMatrix<STATE> > MassMatrix(TPZPoroElasticMF2d * mymaterial, TPZCompMesh* mphysics, int nthreads);
    
    void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
    
    void SetParameters(TPZSimulationData *simulation_data, REAL mod_young, REAL mod_poisson, REAL coef_alpha, REAL coef_Se, REAL permeabil_fluido, REAL visc_fluido, REAL fx, REAL fy,REAL sign);
    
    /*  Malhas geometricas */
    TPZGeoMesh *CreateGMesh();
    
    /*  Uniform refinement */
    void UniformRef(TPZGeoMesh * gmesh, int n_div);

    /* Malhas computacionais */
    TPZCompMesh *CMesh_E(TPZGeoMesh *gmesh, int pOrder, TPZSimulationData *simulation_data); // Malha computacional de elasticidade
    TPZCompMesh *CMesh_q(TPZGeoMesh *gmesh, int pOrder); // Malha computacional de fluxo
    TPZCompMesh *CMesh_p(TPZGeoMesh *gmesh, int pOrder); // Malha computacional de pressão
    TPZCompMesh *CMesh_M(TPZManVector<TPZCompMesh * , 2 > &mesh_vector, TPZGeoMesh *gmesh, int pOrder,TPZSimulationData *simulation_data); // Malha computacional multifísica
    
    //solucao exata
    static void Sol_exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);
    
    //Fractures structure
    void Plot_over_fractures(TPZCompMesh *cmesh);
    void BreakConnectivity(TPZCompMesh &cmesh, std::vector<int> fracture_ids);
    void BreakH1Connectivity(TPZCompMesh &cmesh, std::vector<int> fracture_ids);
    
    //Multiphysics Interfaces
    void AddMultiphysicsInterfaces(TPZCompMesh &cmesh,int mat_frac);
    
    bool insert_fractures_Q = true;
    
    void AdjustIntegrationOrder(TPZSimulationData * sim_data, TPZCompMesh * cmesh_geomechanic, TPZCompMesh * cmesh_reservoir);
    
};


#endif 

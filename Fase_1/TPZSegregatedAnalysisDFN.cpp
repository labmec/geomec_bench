//
//  TPZSegregatedAnalysisDFN.cpp
//  PMRS
//
//  Created by Omar Durán on 9/13/18.
//

#include "TPZSegregatedAnalysisDFN.h"


TPZSegregatedAnalysisDFN::TPZSegregatedAnalysisDFN(){
    m_simulation_data       = NULL;
    m_elastoplast_analysis  = NULL;
    m_darcy_analysis    = NULL;
}

TPZSegregatedAnalysisDFN::~TPZSegregatedAnalysisDFN(){
    
}

TPZSegregatedAnalysisDFN::TPZSegregatedAnalysisDFN(const TPZSegregatedAnalysisDFN & other){
    m_simulation_data       = other.m_simulation_data;
    m_elastoplast_analysis  = other.m_elastoplast_analysis;
    m_darcy_analysis    = other.m_darcy_analysis;
}

void TPZSegregatedAnalysisDFN::ApplyFracMemoryLink(int frac_matid){

    TPZMaterial * frac_material_elastoplast = m_elastoplast_analysis->Mesh()->FindMaterial(frac_matid);
    TPZMaterial * frac_material_darcy = m_darcy_analysis->Mesh()->FindMaterial(frac_matid);
    
    if (!frac_material_elastoplast || !frac_material_darcy) {
        DebugStop();
    }

    TPZMatWithMem<TPZMemoryFracDFN> * matfrac_with_memory_elastoplast = dynamic_cast<TPZMatWithMem<TPZMemoryFracDFN> * >(frac_material_elastoplast);
    TPZMatWithMem<TPZMemoryFracDFN> * matfrac_with_memory_darcy = dynamic_cast<TPZMatWithMem<TPZMemoryFracDFN> * >(frac_material_darcy);
    
    int memdarcy = matfrac_with_memory_darcy->GetMemory()->NElements();
    int memelast = matfrac_with_memory_elastoplast->GetMemory()->NElements();
    if(memdarcy != memelast)
    {
        std::cout << "The integration rules of both meshes are different - bailing out\n";
        DebugStop();
    }
    if (!matfrac_with_memory_elastoplast || !matfrac_with_memory_darcy) {
        DebugStop();
    }
    
    matfrac_with_memory_darcy->SetMemory(matfrac_with_memory_elastoplast->GetMemory());
    
}

void TPZSegregatedAnalysisDFN::ApplyMemoryLink(int matid){
    
    matid = m_simulation_data->Get_elasticity_matid(); //Material Volumétrico

    TPZMaterial * material_elastoplast = m_elastoplast_analysis->Mesh()->FindMaterial(matid);
    TPZMaterial * material_darcy = m_darcy_analysis->Mesh()->FindMaterial(matid);
    
    if (!material_elastoplast || !material_darcy) {
        DebugStop();
    }
    
    TPZMatWithMem<TPZMemoryDFN> * mat_with_memory_elastoplast = dynamic_cast<TPZMatWithMem<TPZMemoryDFN> * >(material_elastoplast);
    TPZMatWithMem<TPZMemoryDFN> * mat_with_memory_darcy = dynamic_cast<TPZMatWithMem<TPZMemoryDFN> * >(material_darcy);
    
    
    if(mat_with_memory_darcy->GetMemory()->NElements() != mat_with_memory_elastoplast->GetMemory()->NElements())
    {
        std::cout << "The integration rules of both meshes are different - bailing out\n";
        DebugStop();
    }
    if (!mat_with_memory_elastoplast || !mat_with_memory_darcy) {
        DebugStop();
    }
    mat_with_memory_darcy->SetMemory(mat_with_memory_elastoplast->GetMemory());
    
}

void TPZSegregatedAnalysisDFN::ConfigurateAnalysis(DecomposeType decompose_E, DecomposeType decompose_M, TPZSimulationData * simulation_data, TPZCompMesh * cmesh_E, TPZCompMesh * cmesh_M, TPZManVector<TPZCompMesh * , 2> & mesh_vec, TPZStack<std::string> & post_pro_var_E, TPZStack<std::string> & post_pro_var_M){
    
    if (!cmesh_E || !cmesh_M) {
        DebugStop();
    }
    
    this->SetSimulationData(simulation_data);
    // totototototo
    bool mustOptimizeBandwidth = false;
    
    m_simulation_data = simulation_data;
    
    // The Geomechanics Simulator
    m_elastoplast_analysis = new TPZPoroElastoPlasticAnalysis;
    m_elastoplast_analysis->SetCompMesh(cmesh_E,mustOptimizeBandwidth);
    m_elastoplast_analysis->ConfigurateAnalysis(decompose_E, m_simulation_data);
    
    // The Reservoir Simulator
    m_darcy_analysis = new TPZDarcyAnalysis;
    m_darcy_analysis->SetCompMesh(cmesh_M,mustOptimizeBandwidth);
    m_darcy_analysis->ConfigurateAnalysis(decompose_M, mesh_vec, m_simulation_data, post_pro_var_M);
    
    for (int imat = 0; imat < simulation_data->Get_volumetric_material_id().size(); imat++) {
        int matid = simulation_data->Get_volumetric_material_id()[imat];
        this->ApplyMemoryLink(matid);
    }
    
    for (int imat_frac = 0; imat_frac < simulation_data->Get_fracture_material_id().size(); imat_frac++) {
        int frac_matid = simulation_data->Get_fracture_material_id()[imat_frac];
        this->ApplyFracMemoryLink(frac_matid);
    }
    
    //this->AdjustIntegrationOrder(cmesh_E,cmesh_M);
    //this->AdjustFractureIntegrationOrder(cmesh_E,cmesh_M);
    
    {
        std::ofstream filecE("MalhaC_E_AfterAdjust.txt"); //Impressão da malha computacional da velocidade (formato txt)
        cmesh_E->Print(filecE);
    }
    {
        std::ofstream filecM("MalhaC_M_AfterAdjust.txt"); //Impressão da malha computacional da velocidade (formato txt)
        cmesh_M->Print(filecM);
    }
    
}


void TPZSegregatedAnalysisDFN::ExecuteOneTimeStep(){
    m_darcy_analysis->ExecuteOneTimeStep();
    if(!m_simulation_data->IsMonoPhasicQ()){
        m_elastoplast_analysis->ExecuteOneTimeStep();
    }
}

void TPZSegregatedAnalysisDFN::PostProcessTimeStep(std::string & geo_file, std::string & res_file){
    m_darcy_analysis->PostProcessTimeStep(res_file);
    if(!m_simulation_data->IsMonoPhasicQ()){
        m_elastoplast_analysis->PostProcessTimeStep(geo_file);
    }
}

void TPZSegregatedAnalysisDFN::ExecuteTimeEvolution(){
//    DebugStop();
    std::string file_darcy("DarcyFlow.vtk");
    std::string file_elastoplast("Elastoplasticity.vtk");

    //Testes
    std::string file_darcy_test("DarcyFlow_test.vtk");
    std::string file_elastoplast_test("Elastoplasticity_teste.vtk");
    
    int n_max_fss_iterations = 10; // @TODO:: MS, please to xml file structure
    int n_enforced_fss_iterations = 1; // @TODO:: MS, please to xml file structure
    int n_time_steps = 1;
    REAL r_norm = m_simulation_data->Get_epsilon_res();
    REAL dx_norm = m_simulation_data->Get_epsilon_cor();
    bool error_stop_criterion_Q = false;
    bool dx_stop_criterion_Q = false;
    this->SetInitialParameters();
    
    for (int it = 0; it < n_time_steps; it++) { //??
        for (int k = 1; k <= n_max_fss_iterations; k++) {
            this->ExecuteOneTimeStep();
            
            if(!m_simulation_data->IsMonoPhasicQ()){
                if (k==1) {
                    this->SetInitialStress();
                }
                this->UpdateParameters();
            }
            error_stop_criterion_Q = (m_darcy_analysis->Get_error() < r_norm) && (m_elastoplast_analysis->Get_error() < r_norm);
            dx_stop_criterion_Q = (m_darcy_analysis->Get_dx_norm() < dx_norm) && (m_elastoplast_analysis->Get_dx_norm() < dx_norm);
            this->PostProcessTimeStep(file_elastoplast_test, file_darcy_test);
            
            if ((error_stop_criterion_Q && (k > n_enforced_fss_iterations)) || dx_stop_criterion_Q) {
                this->PostProcessTimeStep(file_elastoplast, file_darcy);
                std::cout << "TPZSegregatedAnalysisDFN:: Iterative process converged with residue norm for Darcy = " << m_darcy_analysis->Get_error() << std::endl;
                std::cout << "TPZSegregatedAnalysisDFN:: Iterative process converged with residue norm for Elastoplasticity = " << m_elastoplast_analysis->Get_error() << std::endl;
                UpdateState();
                break;
            }
        }
    }

}

void TPZSegregatedAnalysisDFN::SetInitialStress(){
    
    // Updating volumetric parameters :
    
    int matid = m_simulation_data->Get_elasticity_matid();
    
    TPZMaterial * material_elastoplast = m_elastoplast_analysis->Mesh()->FindMaterial(matid);
    
    TPZMatWithMem<TPZMemoryDFN> * mat_with_memory_elastoplast = dynamic_cast<TPZMatWithMem<TPZMemoryDFN> * >(material_elastoplast);
    
    long N_ipoints = mat_with_memory_elastoplast->GetMemory().get()->NElements();
    for (int ip_index = 0 ; ip_index < N_ipoints; ip_index++) {
        TPZMemoryDFN & memory = mat_with_memory_elastoplast->GetMemory().get()->operator[](ip_index);
        memory.SetSigma_0(memory.GetSigma_n());
    }
}


void TPZSegregatedAnalysisDFN::SetInitialParameters(){
    
    // Updating volumetric parameters :
    
    int matid = m_simulation_data->Get_elasticity_matid();
    
    TPZMaterial * material_elastoplast = m_elastoplast_analysis->Mesh()->FindMaterial(matid);
    TPZMaterial * material_darcy = m_darcy_analysis->Mesh()->FindMaterial(matid);
    
    if (!material_elastoplast || !material_darcy) {
        DebugStop();
    }
    
    TPZMatWithMem<TPZMemoryDFN> * mat_with_memory_darcy = dynamic_cast<TPZMatWithMem<TPZMemoryDFN> * >(material_darcy);
    
    long N_ipoints = mat_with_memory_darcy->GetMemory().get()->NElements();
    for (int ip_index = 0 ; ip_index < N_ipoints; ip_index++) {
        TPZMemoryDFN & memory = mat_with_memory_darcy->GetMemory().get()->operator[](ip_index);
        TPZFNMatrix<9,REAL>  k_0 = m_simulation_data->Get_PermeabilityTensor_0();
        REAL phi0 = m_simulation_data->Get_Porosity_0();
        memory.Setkappa_0(k_0);
        memory.Setkappa(k_0);
        memory.Setphi_0(phi0);
    }
    
    // Updating Fracture parameters :
    
    for (int imat_frac = 0; imat_frac < m_simulation_data->Get_fracture_material_id().size(); imat_frac++) {
        int frac_matid = m_simulation_data->Get_fracture_material_id()[imat_frac];
        
        TPZMaterial * frac_material_elastoplast = m_elastoplast_analysis->Mesh()->FindMaterial(frac_matid);
        TPZMaterial * frac_material_darcy = m_darcy_analysis->Mesh()->FindMaterial(frac_matid);
        
        
        
        if (!frac_material_elastoplast || !frac_material_darcy) {
            DebugStop();
        }
        
        TPZMatWithMem<TPZMemoryFracDFN> * matfrac_with_memory_darcy = dynamic_cast<TPZMatWithMem<TPZMemoryFracDFN> * >(frac_material_darcy);
        
        long N_ipoints = matfrac_with_memory_darcy->GetMemory().get()->NElements();
        
        for (int ip_index = 0 ; ip_index < N_ipoints; ip_index++) {
            TPZMemoryFracDFN & memory_frac = matfrac_with_memory_darcy->GetMemory().get()->operator[](ip_index);
            TPZFNMatrix<9,REAL>  k_0(3,3,0.);
            k_0(0,0) = m_simulation_data->Get_Permeability_0();
            
            REAL Vm = m_simulation_data->Get_Vm().find(frac_matid)->second; //Max fracture closure
            REAL a0 = m_simulation_data->Get_a0().find(frac_matid)->second;
            REAL Kni = m_simulation_data->Get_Kni().find(frac_matid)->second;
            
            REAL Du_0 = Vm-a0; // initial closure
            
            Du_0 = (-30.*Vm)/(-30.-Kni*Vm);
            
            memory_frac.SetVm(Vm);
            memory_frac.SetDu_0(Du_0);
            memory_frac.Setkappa_0(k_0);
            memory_frac.Setkappa(k_0);
        }
        
        
    }
    
}



void TPZSegregatedAnalysisDFN::UpdateParameters(){
    
    // Updating volumetric parameters :
    
    int matid = m_simulation_data->Get_elasticity_matid();
    
    TPZMaterial * material_elastoplast = m_elastoplast_analysis->Mesh()->FindMaterial(matid);
    TPZMaterial * material_darcy = m_darcy_analysis->Mesh()->FindMaterial(matid);
    
    if (!material_elastoplast || !material_darcy) {
        DebugStop();
    }
    
    TPZMatWithMem<TPZMemoryDFN> * mat_with_memory_elastoplast = dynamic_cast<TPZMatWithMem<TPZMemoryDFN> * >(material_elastoplast);
    TPZMatWithMem<TPZMemoryDFN> * mat_with_memory_darcy = dynamic_cast<TPZMatWithMem<TPZMemoryDFN> * >(material_darcy);
    
    long N_ipoints = mat_with_memory_darcy->GetMemory().get()->NElements();
    
    for (int ip_index = 0 ; ip_index < N_ipoints; ip_index++) {
        TPZMemoryDFN & memory = mat_with_memory_darcy->GetMemory().get()->operator[](ip_index);
        TPZFNMatrix<9,REAL>  k_0 = memory.kappa_0();
        REAL phi_0 = memory.phi_0();
        REAL nu = m_simulation_data->Get_Poisson();
        REAL E = m_simulation_data->Get_Eyoung();
        TPZFNMatrix<9,REAL> k_n(3,3,0.);
        for(int i = 0; i < k_0.Rows(); i++){
            k_n(i,i) = memory.Permeability(k_0(i,i),phi_0,nu,E);
        }
        memory.Setkappa(k_n);
    }
    
    // Updating Fracture parameters :
    
    for (int imat_frac = 0; imat_frac < m_simulation_data->Get_fracture_material_id().size(); imat_frac++) {
        int frac_matid = m_simulation_data->Get_fracture_material_id()[imat_frac];
        
        TPZMaterial * frac_material_elastoplast = m_elastoplast_analysis->Mesh()->FindMaterial(frac_matid);
        TPZMaterial * frac_material_darcy = m_darcy_analysis->Mesh()->FindMaterial(frac_matid);
        
        if (!frac_material_elastoplast || !frac_material_darcy) {
            DebugStop();
        }
        
        TPZMatWithMem<TPZMemoryFracDFN> * matfrac_with_memory_elastoplast = dynamic_cast<TPZMatWithMem<TPZMemoryFracDFN> * >(frac_material_elastoplast);
        TPZMatWithMem<TPZMemoryFracDFN> * matfrac_with_memory_darcy = dynamic_cast<TPZMatWithMem<TPZMemoryFracDFN> * >(frac_material_darcy);

        long N_ipoints = matfrac_with_memory_darcy->GetMemory().get()->NElements();
        
        for (int ip_index = 0 ; ip_index < N_ipoints; ip_index++) {
            TPZMemoryFracDFN & memory_frac = matfrac_with_memory_darcy->GetMemory().get()->operator[](ip_index);
            TPZFNMatrix<9,REAL>  k_0 = memory_frac.kappa_0();
            
            TPZFNMatrix<9,REAL> k_n(k_0.Rows(),k_0.Cols(),0.);
            for(int i = 0; i < k_0.Rows(); i++){
                k_n(i,i) = memory_frac.Permeability(k_0(i,i));
            }
            memory_frac.Setkappa(k_n);
        }
        

    }
    
}




void TPZSegregatedAnalysisDFN::UpdateState(){
    m_darcy_analysis->UpdateState();
    m_elastoplast_analysis->UpdateState();
}

void TPZSegregatedAnalysisDFN::AdjustIntegrationOrder(TPZCompMesh * cmesh_o, TPZCompMesh * cmesh_d){
    
    // Assuming the cmesh_o as directive.
    
    cmesh_d->LoadReferences();
    int nel_o = cmesh_o->NElements();
    int nel_d = cmesh_d->NElements();
    
    int nvol_el_o = 0;
    int nvol_el_d = 0;
    
    for (long el = 0; el < nel_o; el++) {
        TPZCompEl *cel_o = cmesh_o->Element(el);
        int matid_o = cel_o->Material()->Id();
        for (int imat = 0; imat < m_simulation_data->Get_volumetric_material_id().size(); imat++) {
            int matid = m_simulation_data->Get_volumetric_material_id()[imat];
            if (matid==matid_o) {
                nvol_el_o ++;
                break;
            }
        }
    }
    
    for (long el = 0; el < nel_d; el++) {
        TPZCompEl *cel_d = cmesh_d->Element(el);
        int matid_d = cel_d->Material()->Id();
        for (int imat = 0; imat < m_simulation_data->Get_volumetric_material_id().size(); imat++) {
            int matid = m_simulation_data->Get_volumetric_material_id()[imat];
            if (matid==matid_d) {
                nvol_el_d ++;
                break;
            }
        }
    }
    
    if (nvol_el_o != nvol_el_d) {
        std::cout << "The geometrical partitions are not the same." << std::endl;
        DebugStop();
    }
    
    for (long el = 0; el < nel_o; el++) {
        TPZCompEl *cel_o = cmesh_o->Element(el);
        if (!cel_o) {
            continue;
        }
        
        TPZGeoEl * gel = cel_o->Reference();
        if (!gel) {
            continue;
        }
        
        // Finding the other computational element
        TPZCompEl * cel_d = gel->Reference();
        if (!cel_d) {
            continue;
        }
        cel_o->SetFreeIntPtIndices();
        cel_o->ForcePrepareIntPtIndices();
        const TPZIntPoints & rule = cel_o->GetIntegrationRule();
        TPZIntPoints * cloned_rule = rule.Clone();
        
        TPZManVector<int64_t,20> indices;
        cel_o->GetMemoryIndices(indices);
        cel_d->SetFreeIntPtIndices();
        cel_d->SetMemoryIndices(indices);
        cel_d->SetIntegrationRule(cloned_rule);
    }
    
#ifdef PZDEBUG
    std::ofstream out_geo("Cmesh_origin_adjusted.txt");
    cmesh_o->Print(out_geo);
#endif
    
    
#ifdef PZDEBUG
    std::ofstream out_res("Cmesh_destination_adjusted.txt");
    cmesh_d->Print(out_res);
#endif
    
}

void TPZSegregatedAnalysisDFN::AdjustFractureIntegrationOrder(TPZCompMesh * cmesh_o, TPZCompMesh * cmesh_d){
    
    // Assuming the cmesh_o as directive.
    
    cmesh_d->LoadReferences();
    int nel_o = cmesh_o->NElements();
    int nel_d = cmesh_d->NElements();
    
    int nvol_el_o = 0;
    int nvol_el_d = 0;
    
    for (long el = 0; el < nel_o; el++) {
        TPZCompEl *cel_o = cmesh_o->Element(el);
        int matid_o = cel_o->Material()->Id();
        for (int imat = 0; imat < m_simulation_data->Get_fracture_material_id().size(); imat++) {
            int matid = m_simulation_data->Get_fracture_material_id()[imat];
            if (matid==matid_o) {
                nvol_el_o ++;
                break;
            }
        }
    }
    
    for (long el = 0; el < nel_d; el++) {
        TPZCompEl *cel_d = cmesh_d->Element(el);
        int matid_d = cel_d->Material()->Id();
        for (int imat = 0; imat < m_simulation_data->Get_fracture_material_id().size(); imat++) {
            int matid = m_simulation_data->Get_fracture_material_id()[imat];
            if (matid==matid_d) {
                nvol_el_d ++;
                break;
            }
        }
    }
    
    if (nvol_el_o != nvol_el_d*2) {
        std::cout << "The geometrical partitions are not the same." << std::endl;
        DebugStop();
    }
    
    for (long el = 0; el < nel_o; el++) {
        TPZCompEl *cel_o = cmesh_o->Element(el);
        if (!cel_o) {
            continue;
        }
        
        TPZGeoEl * gel = cel_o->Reference();
        if (!gel) {
            continue;
        }
        
        // Finding the other computational element
        TPZCompEl * cel_d = gel->Reference();
        if (!cel_d) {
            continue;
        }
        cel_o->SetFreeIntPtIndices();
        cel_o->ForcePrepareIntPtIndices();
        const TPZIntPoints & rule = cel_o->GetIntegrationRule();
        TPZIntPoints * cloned_rule = rule.Clone();
        
        TPZManVector<int64_t,20> indices;
        cel_o->GetMemoryIndices(indices);
        cel_d->SetFreeIntPtIndices();
        cel_d->SetMemoryIndices(indices);
        cel_d->SetIntegrationRule(cloned_rule);
    }
    
#ifdef PZDEBUG
    std::ofstream out_geo("Cmesh_origin_adjusted.txt");
    cmesh_o->Print(out_geo);
#endif
    
    
#ifdef PZDEBUG
    std::ofstream out_res("Cmesh_destination_adjusted.txt");
    cmesh_d->Print(out_res);
#endif
    
}

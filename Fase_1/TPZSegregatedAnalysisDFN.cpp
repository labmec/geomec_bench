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

void TPZSegregatedAnalysisDFN::ApplyMemoryLink(){
    
    int matid = 1; //Material Volumétrico
        
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
    bool mustOptimizeBandwidth = true;
    
    // The Geomechanics Simulator
    m_elastoplast_analysis = new TPZMatElastoPlasticAnalysis;
    m_elastoplast_analysis->SetCompMesh(cmesh_E,mustOptimizeBandwidth);
    m_elastoplast_analysis->ConfigurateAnalysis(decompose_E, m_simulation_data);
    
    // The Reservoir Simulator
    m_darcy_analysis = new TPZDarcyAnalysis;
    m_darcy_analysis->SetCompMesh(cmesh_M,mustOptimizeBandwidth);
    m_darcy_analysis->ConfigurateAnalysis(decompose_M, mesh_vec, m_simulation_data, post_pro_var_M);
    
    this->ApplyMemoryLink();
    
}


void TPZSegregatedAnalysisDFN::ExecuteOneTimeStep(){
    m_darcy_analysis->ExecuteOneTimeStep();
    m_elastoplast_analysis->ExecuteOneTimeStep();
}

void TPZSegregatedAnalysisDFN::PostProcessTimeStep(std::string & geo_file, std::string & res_file){
    m_darcy_analysis->PostProcessTimeStep(res_file);
    m_elastoplast_analysis->PostProcessTimeStep(geo_file);
}

void TPZSegregatedAnalysisDFN::ExecuteTimeEvolution(){
//    DebugStop();
    std::string file_darcy("DarcyFlow.vtk");
    std::string file_elastoplast("Elastoplasticity.vtk");

    int n_max_fss_iterations = 10; // @TODO:: MS, please to xml file structure
    int n_enforced_fss_iterations = 1; // @TODO:: MS, please to xml file structure
    int n_time_steps = 1;
    REAL r_norm = m_simulation_data->Get_epsilon_res();
    REAL dx_norm = m_simulation_data->Get_epsilon_cor();
    bool error_stop_criterion_Q = false;
    bool dx_stop_criterion_Q = false;
    for (int it = 0; it < n_time_steps; it++) {
        for (int k = 1; k <= n_max_fss_iterations; k++) {
            this->ExecuteOneTimeStep();
            error_stop_criterion_Q = (m_darcy_analysis->Get_error() < r_norm) && (m_elastoplast_analysis->Get_error() < r_norm);
            dx_stop_criterion_Q = (m_darcy_analysis->Get_dx_norm() < dx_norm) && (m_elastoplast_analysis->Get_dx_norm() < dx_norm);

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

void TPZSegregatedAnalysisDFN::UpdateState(){
    m_darcy_analysis->UpdateState();
    m_elastoplast_analysis->UpdateState();
}

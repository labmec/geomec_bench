//
//  TPZDarcyAnalysis.cpp
//  Benchmark0a
//
//  Created by Pablo Carvalho on 14/09/18.
//

#include "TPZDarcyAnalysis.h"
#include "pzbndcond.h"
#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("Benchmark.DarcyAnalysis"));
#endif
TPZDarcyAnalysis::TPZDarcyAnalysis() : TPZAnalysis(){
    
    m_simulation_data = NULL;
    m_X_n.Resize(0, 0);
    m_X.Resize(0, 0);
    m_mesh_vec.Resize(0);
    m_error = 0;
    m_dx_norm = 0;
    m_k_iterations = 0;
    m_post_processor = NULL;
    m_var_names.resize(0);
    m_vec_var_names.resize(0);
    
}

TPZDarcyAnalysis::~TPZDarcyAnalysis(){
    
}

TPZDarcyAnalysis::TPZDarcyAnalysis(const TPZDarcyAnalysis & other){
    
    m_simulation_data   = other.m_simulation_data;
    m_X_n               = other.m_X_n;
    m_X                 = other.m_X;
    m_mesh_vec          = other.m_mesh_vec;
    m_error             = other.m_error;
    m_dx_norm           = other.m_dx_norm;
    m_k_iterations      = other.m_k_iterations;
    m_post_processor    = other.m_post_processor;
    m_var_names         = other.m_var_names;
    m_vec_var_names     = other.m_vec_var_names;
    
}

void TPZDarcyAnalysis::ConfigurateAnalysis(DecomposeType decomposition, TPZManVector<TPZCompMesh * , 2> & mesh_vec, TPZSimulationData * simulation_data, TPZVec<std::string> & var_names){
    
    SetSimulationData(simulation_data);
    TPZStepSolver<STATE> step;
    unsigned int n_threads = m_simulation_data->Get_n_threads();

    if(!Mesh()){
        std::cout << "Call SetCompMesh method." << std::endl;
        DebugStop();
    }
    
    m_mesh_vec = mesh_vec;
    switch (decomposition) {
        case ELU:
        {
//#ifdef USING_MKL
//            TPZSpStructMatrix struct_mat(Mesh());
//            struct_mat.SetNumThreads(n_threads);
//            this->SetStructuralMatrix(struct_mat);
//#else
            TPZSkylineNSymStructMatrix struct_mat(Mesh());
            struct_mat.SetNumThreads(n_threads);
            this->SetStructuralMatrix(struct_mat);
//#endif
        }
            break;
        case ELDLt:
        {
            TPZSymetricSpStructMatrix struct_mat(Mesh());
            struct_mat.SetNumThreads(n_threads);
            this->SetStructuralMatrix(struct_mat);
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
    step.SetDirect(decomposition);
    this->SetSolver(step);
    this->Solution().Resize(Mesh()->Solution().Rows(), 1);
    m_X.Resize(Mesh()->Solution().Rows(), 1);
    m_X_n.Resize(Mesh()->Solution().Rows(), 1);

    m_post_processor = new TPZPostProcAnalysis;
    m_post_processor->SetCompMesh(Mesh());

    int n_vols = m_simulation_data->Get_volumetric_material_id().size();
    TPZManVector<int,10> post_mat_id(n_vols);
    for (int ivol = 0; ivol < n_vols; ivol++)
    {
        int matid = m_simulation_data->Get_volumetric_material_id()[ivol];
        post_mat_id[ivol] = matid;
    }

    for (auto i : var_names) {
        m_var_names.Push(i);
    }

    m_post_processor->SetPostProcessVariables(post_mat_id, m_var_names);
    int dim = Mesh()->Dimension();
    int div = 0;
    TPZStack< std::string> vecnames;
    std::string plotfile("Benchmark_Mono_DarcyTest.vtk");

    m_post_processor->DefineGraphMesh(dim,m_var_names,vecnames,plotfile);

    TPZFStructMatrix structmatrix(m_post_processor->Mesh());
    structmatrix.SetNumThreads(n_threads);
    m_post_processor->SetStructuralMatrix(structmatrix);

}

void TPZDarcyAnalysis::ExecuteNewtonInteration(){
    this->Assemble();
    this->Rhs() *= -1.0;
    this->Solve();
}

void TPZDarcyAnalysis::ExecuteOneTimeStep(){
  
    // m_X means the solution at the previous time step
    if (m_simulation_data->IsInitialStateQ()) {
        m_X = Solution();
    }
    // Set current state false means overwriting p of the memory
    m_simulation_data->SetCurrentStateQ(false);
    // Accect time solution means writing one of the vectors of this object in the memory
    AcceptTimeStepSolution();
    
    //    // Initial guess
    //    m_X_n = m_X;
    // Set current state true means overwriting p_n of the memory object
    m_simulation_data->SetCurrentStateQ(true);
    
    // Accept time solution here means writing one of the vectors of the object into the memory
    this->AcceptTimeStepSolution();
    
    
    std::ofstream plotDarcyEK("DarcyStiffness.txt");
    std::ofstream plotDarcyEF("DarcyRhs.txt");
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        fRhs.Print("Rhs =",sout);
        //EFormatted, EInputFormat, EMathematicaInput, EMatlabNonZeros, EMatrixMarket
      //  fSolver->Matrix()->Print("ek = ",plotDarcyEK,EMathematicaInput);
        fRhs.Print("ef = ",plotDarcyEF,EMathematicaInput);
        
        PrintVectorByElement(sout, fRhs);
        PrintVectorByElement(sout, fSolution);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    TPZFMatrix<STATE> dx;
    bool residual_stop_criterion_Q = false;
    bool correction_stop_criterion_Q = false;
    REAL norm_res, norm_dx;
    REAL r_norm = m_simulation_data->Get_epsilon_res();
    REAL dx_norm = m_simulation_data->Get_epsilon_cor();
    int n_it = m_simulation_data->Get_n_iterations();
    
    for (int i = 1; i <= n_it; i++) {
        this->ExecuteNewtonInteration();
        
        dx = Solution();
        norm_dx  = Norm(dx);
        m_X_n += dx;
        LoadCurrentState();
        
        AssembleResidual();

        // Test: Changing boundary condition into Dirichlet
//        std::map<int, TPZMaterial * >::iterator matit;
//        for(matit = fCompMesh->MaterialVec().begin(); matit != fCompMesh->MaterialVec().end(); matit++)
//        {
//
//            TPZMaterial *mat = matit->second;
//            if(mat->Id()==503){
//                TPZBndCond *bc = dynamic_cast<TPZBndCond *> (mat);
//                bc->SetType(0);
//            }
//        }
//        AssembleResidual();

        
#ifdef LOG4CXX
       if(logger->isDebugEnabled())
       {
           std::stringstream sout;
           fRhs.Print("Rhs =",sout);
           {
           std::ofstream plotDarcyEK("DarcyStiffness.txt");
           std::ofstream plotDarcyEF("DarcyRhs.txt");
           fSolver->Matrix()->Print("ek = ",plotDarcyEK,EMathematicaInput);
           fRhs.Print("ef = ",plotDarcyEF,EMathematicaInput);
           }
           PrintVectorByElement(sout, fRhs);
           PrintVectorByElement(sout, fSolution);
           LOGPZ_DEBUG(logger, sout.str())
       }
#endif
        
        norm_res = Norm(Rhs());
        residual_stop_criterion_Q   = norm_res < r_norm;
        correction_stop_criterion_Q = norm_dx  < dx_norm;
        
        m_k_iterations = i;
        m_error = norm_res;
        m_dx_norm = norm_dx;
        
        

        if (residual_stop_criterion_Q &&  correction_stop_criterion_Q) {
#ifdef PZDEBUG
            std::cout << "TPMRSMonoPhasicAnalysis:: Nonlinear process converged with residue norm = " << norm_res << std::endl;
            std::cout << "TPMRSMonoPhasicAnalysis:: Number of iterations = " << i << std::endl;
            std::cout << "TPMRSMonoPhasicAnalysis:: Correction norm = " << norm_dx << std::endl;
#endif
            this->AcceptTimeStepSolution();
            break;
        }
    }
    
    if (residual_stop_criterion_Q == false) {
        std::cout << "TPMRSMonoPhasicAnalysis:: Nonlinear process not converged with residue norm = " << norm_res << std::endl;
    }
    
  
}

void TPZDarcyAnalysis::UpdateState(){
    m_X = m_X_n;
}

void TPZDarcyAnalysis::PostProcessTimeStep(std::string & file, bool is_stantdard_post_pro_Q){
    
    if (is_stantdard_post_pro_Q) {
        this->StandardPostProcessTimeStep(file);
        return;
    }
    
    int dim = Mesh()->Dimension();
    int div = 0;
    TPZStack< std::string> vecnames;
    m_post_processor->TransferSolution();
    m_post_processor->DefineGraphMesh(dim,m_var_names,vecnames,file);
    m_post_processor->PostProcess(div,dim);
}

void TPZDarcyAnalysis::StandardPostProcessTimeStep(std::string & file){
    
    int postProcessResolution = m_simulation_data->Get_vtk_resolution();
    int dim = Mesh()->Reference()->Dimension();
    this->DefineGraphMesh(dim,m_var_names,m_vec_var_names,file);
    this->PostProcess(postProcessResolution,dim);
    
}

void TPZDarcyAnalysis::AcceptTimeStepSolution(){

  //  DebugStop();
    
    bool state = m_simulation_data->IsCurrentStateQ();
    if (state) {
        
        // must accept solution changes a global data structure shared by the material objects
        // which indicates the solution should be overwritten in memory
        m_simulation_data->Set_must_accept_solution_Q(true);
        // load current state copies m_X_n into the solution vector
        LoadCurrentState();
        // puts the solution vector into a variable depending on yet another global variable
        AssembleResidual();
        m_simulation_data->Set_must_accept_solution_Q(false);
    }else{
        // m_simulation_data is pointer shared by the material object
        // this call forces the solution to be loaded into the memory object
        m_simulation_data->Set_must_accept_solution_Q(true);
        // put m_X in the mesh solution
        LoadLastState();
        // copy the state vector into the memory because must_accept_solution_Q in the m_simulation_data is true
        AssembleResidual();
        m_simulation_data->Set_must_accept_solution_Q(false);
    }
}


void TPZDarcyAnalysis::LoadCurrentState(){
    LoadSolution(m_X_n);
    if(m_simulation_data->Get_is_dual_formulation_Q()){
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_mesh_vec, Mesh());
    }
}

void TPZDarcyAnalysis::LoadLastState(){
    LoadSolution(m_X);
    if(m_simulation_data->Get_is_dual_formulation_Q()){
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_mesh_vec, Mesh());
    }
}


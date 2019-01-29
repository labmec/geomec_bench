//
//  TPZPoroElastoPlasticAnalysis.hpp
//  Benchmark0a
//
//  Created by Pablo Carvalho on 14/09/18.
//


#include "TPZPoroElastoPlasticAnalysis.h"
#include "pzbndcond.h"
#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr loggerElast(Logger::getLogger("Benchmark.ElastAnalysis"));
#endif


TPZPoroElastoPlasticAnalysis::TPZPoroElastoPlasticAnalysis() : TPZAnalysis(){
    
    m_simulation_data = NULL;
    m_X_n.Resize(0, 0);
    m_X.Resize(0, 0);
    m_error = 0;
    m_dx_norm = 0;
    m_k_iterations = 0;
    m_post_processor = NULL;
    m_var_names.resize(0);
    
}

TPZPoroElastoPlasticAnalysis::~TPZPoroElastoPlasticAnalysis(){
    
}

TPZPoroElastoPlasticAnalysis::TPZPoroElastoPlasticAnalysis(const TPZPoroElastoPlasticAnalysis & other){
    
    m_simulation_data   = other.m_simulation_data;
    m_X_n               = other.m_X_n;
    m_X                 = other.m_X;
    m_error             = other.m_error;
    m_dx_norm           = other.m_dx_norm;
    m_k_iterations      = other.m_k_iterations;
    m_post_processor    = other.m_post_processor;
    m_var_names         = other.m_var_names;
    
}

void TPZPoroElastoPlasticAnalysis::ConfigurateAnalysis(DecomposeType decomposition, TPZSimulationData * simulation_data){
    SetSimulationData(simulation_data);
    TPZStepSolver<STATE> step;
    unsigned int number_threads = m_simulation_data->Get_n_threads();
    
    if(!Mesh()){
        std::cout << "Call SetCompMesh method." << std::endl;
        DebugStop();
    }
    
    switch (decomposition) {
        case ECholesky:
        {
            TPZSkylineStructMatrix struct_mat(Mesh());
            struct_mat.SetNumThreads(number_threads);
            this->SetStructuralMatrix(struct_mat);
        }
            break;
        case ELDLt:
        {
//            TPZSymetricSpStructMatrix struct_mat(Mesh());
            TPZSkylineStructMatrix struct_mat(Mesh());
            struct_mat.SetNumThreads(number_threads);
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
    
    int n_threads = m_simulation_data->Get_n_threads();
    m_post_processor = new TPZPostProcAnalysis;
    m_post_processor->SetCompMesh(Mesh());
    
    TPZManVector<std::pair<int, TPZManVector<int,12>>,12>  material_ids = m_simulation_data->MaterialIds();
    TPZManVector<int,10> post_mat_id(1);
    
    int matid = 1;
    post_mat_id[0] = matid;
    
    // @TODO:: MS, please transfer from xml file
    m_var_names.Push("ux");
    m_var_names.Push("uy");
    m_var_names.Push("sxx");
    m_var_names.Push("syy");
    m_var_names.Push("szz");
    m_var_names.Push("exx");
    m_var_names.Push("eyy");
    m_var_names.Push("ezz");
    m_var_names.Push("epxx");
    m_var_names.Push("epyy");
    m_var_names.Push("epzz");
    
//    if (m_simulation_data->Dimension() == 3) {
//        m_var_names.Push("uz");
//    }
    
    m_post_processor->SetPostProcessVariables(post_mat_id, m_var_names);
    
    TPZFStructMatrix structmatrix(m_post_processor->Mesh());
    structmatrix.SetNumThreads(n_threads);
    m_post_processor->SetStructuralMatrix(structmatrix);
    
}

void TPZPoroElastoPlasticAnalysis::ExecuteNewtonInteration(){
    this->Assemble();
    this->Rhs() *= -1.0;
    {
        std::ofstream out("RhsExecuteNewton.txt");
        PrintVectorByElement(out, this->Rhs());
    }
    this->Solve();
    {
        std::ofstream out("SolExecuteNewton.txt");
        PrintVectorByElement(out, this->Solution());
    }
}

void TPZPoroElastoPlasticAnalysis::ExecuteOneTimeStep(bool must_accept_solution_Q){
    
    if (m_simulation_data->IsInitialStateQ()) {
        m_X = Solution();
    }
    
    m_simulation_data->SetCurrentStateQ(false);
    m_simulation_data->SetInitialStateQ(true);
    AcceptPseudoTimeStepSolution();
    
    m_simulation_data->SetInitialStateQ(false);
    m_simulation_data->SetCurrentStateQ(true);
    AcceptPseudoTimeStepSolution();
    //    // Reset du to zero
    //    Solution().Zero();
    //    LoadSolution(Solution());
    
    TPZFMatrix<STATE> desloc(Solution());
    
    std::ofstream plotElasticEK("ElastStiffness.txt");
    std::ofstream plotElasticEF("ElastRhs.txt");
#ifdef LOG4CXX
    if(loggerElast->isDebugEnabled())
    {
        std::stringstream sout;
        fRhs.Print("Rhs =",sout);
        
      //  fSolver->Matrix()->Print("ek = ",plotElasticEK,EMathematicaInput);
        fRhs.Print("ef = ",plotElasticEF,EMathematicaInput);
        
        PrintVectorByElement(sout, fRhs);
        PrintVectorByElement(sout, fSolution);
        LOGPZ_DEBUG(loggerElast, sout.str())
        m_post_processor->TransferSolution();
        m_post_processor->PostProcess(0);
    }
#endif
    
    
    
    
    bool residual_stop_criterion_Q = false;
    bool correction_stop_criterion_Q = false;
    REAL norm_res, norm_desloc;
    REAL r_norm = m_simulation_data->Get_epsilon_res();
    REAL dx_norm = m_simulation_data->Get_epsilon_cor();
    int n_it = m_simulation_data->Get_n_iterations();
    
    for (int i = 1; i <= n_it; i++) {
        this->ExecuteNewtonInteration();
        desloc = Solution();
        norm_desloc  = Norm(Solution());
        //LoadSolution(desloc);
        m_X_n += desloc;
        LoadCurrentState();
        AssembleResidual();
    
        
        std::ofstream plotElasticEK("ElastStiffness.txt");
        std::ofstream plotElasticEF("ElastRhs.txt");
#ifdef LOG4CXX
        if(loggerElast->isDebugEnabled())
        {
            std::stringstream sout;
            fRhs.Print("Rhs =",sout);
            
            fSolver->Matrix()->Print("ek = ",plotElasticEK,EMathematicaInput);
            fRhs.Print("ef = ",plotElasticEF,EMathematicaInput);
            
            PrintVectorByElement(sout, fRhs);
            PrintVectorByElement(sout, fSolution);
            LOGPZ_DEBUG(loggerElast, sout.str())
            m_post_processor->TransferSolution();
            m_post_processor->PostProcess(0);
        }
#endif
        std::string file_elastoplast_test("Elastoplasticity_teste_before.vtk");
        this->PostProcessTimeStep(file_elastoplast_test);
        
        norm_res = Norm(this->Rhs());
        residual_stop_criterion_Q   = norm_res < r_norm;
        correction_stop_criterion_Q = norm_desloc  < dx_norm;
        
        m_k_iterations = i;
        m_error = norm_res;
        m_dx_norm = norm_desloc;
        
        
        if (residual_stop_criterion_Q && correction_stop_criterion_Q) {
#ifdef PZDEBUG
            std::cout << "TPZPoroElastoPlasticAnalysis:: Nonlinear process converged with residue norm = " << norm_res << std::endl;
            std::cout << "TPZPoroElastoPlasticAnalysis:: Number of iterations = " << i << std::endl;
            std::cout << "TPZPoroElastoPlasticAnalysis:: Correction norm = " << norm_desloc << std::endl;
#endif
            this->AcceptPseudoTimeStepSolution();
#ifdef PZDEBUG
            
            //Solution().Zero();
            //LoadSolution();
            //m_X_n = Solution();
//            LoadCurrentState();
//            //AssembleResidual();
//            AssembleResidual();
//            {
//            std::ofstream fileRhs("SaidaRHS.txt");
//            PrintVectorByElement(fileRhs, this->Rhs());
//            }
//            REAL norm_res_c = Norm(this->Rhs());
//            std::cout << " norm = " << norm_res_c << std::endl;
            

            /// Print Interface Memory vol. elements
            std::ofstream fileMemVol("MemoryDFN.txt", std::ofstream::app);
            fileMemVol <<"For porous media :" << std::endl;
            int mat_ID = m_simulation_data->Get_volumetric_material_id()[0];
            TPZMaterial * mat_volumetric = this->Mesh()->FindMaterial(mat_ID);
            TPZMatWithMem<TPZMemoryDFN> * mat_with_volmem = dynamic_cast<TPZMatWithMem<TPZMemoryDFN> * >(mat_volumetric);
           // mat_with_volmem->GetMemory().Print(fileMemVol);
            mat_with_volmem->Print(fileMemVol);
            
            /// Print Interface Memory Left and Right
            
            if(m_simulation_data->Get_insert_fractures_Q()){
                std::ofstream fileMemInter("MemoryInterfaces.txt", std::ofstream::app);
                fileMemInter <<"For Interfaces left :" << std::endl;
                int interLeft_ID = m_simulation_data->Get_interfaceLeft_id()[0];
                TPZMaterial * mat_interLeft = this->Mesh()->FindMaterial(interLeft_ID);
                TPZMatWithMem<TPZInterfaceMemory, TPZDiscontinuousGalerkin> * mat_with_interLeft = dynamic_cast<TPZMatWithMem<TPZInterfaceMemory, TPZDiscontinuousGalerkin> * >(mat_interLeft);
                mat_with_interLeft->Print(fileMemInter);
                
                fileMemInter <<"For Interfaces right :"<< i << std::endl;
                int interRight_ID = m_simulation_data->Get_interfaceRight_id()[0];
                TPZMaterial * mat_interRight = this->Mesh()->FindMaterial(interRight_ID);
                TPZMatWithMem<TPZInterfaceMemory, TPZDiscontinuousGalerkin> * mat_with_interRight = dynamic_cast<TPZMatWithMem<TPZInterfaceMemory, TPZDiscontinuousGalerkin> * >(mat_interRight);
                mat_with_interRight->Print(fileMemInter);
            
                /// Print Fracture Memory
                std::ofstream fileMemFrac("MemoryFrac.txt", std::ofstream::app);
                int n_fracs = m_simulation_data->Get_fracture_material_id().size();
                for(int i = 0; i< n_fracs; i++) {
                    fileMemFrac <<"For fracture number = "<< i << std::endl;
                    int frac_ID = m_simulation_data->Get_fracture_material_id()[i];
                    TPZMaterial * material_frac = this->Mesh()->FindMaterial(frac_ID);
                    TPZMatWithMem<TPZMemoryFracDFN> * mat_with_memory_frac = dynamic_cast<TPZMatWithMem<TPZMemoryFracDFN> * >(material_frac);
                    mat_with_memory_frac->Print(fileMemFrac);
                }
            }
//            memory.Print(fileMemFrac);
            
#endif
            break;
        }
    }
    
    if (residual_stop_criterion_Q == false) {
        std::cout << "TPZPoroElastoPlasticAnalysis:: Nonlinear process not converged with residue norm = " << norm_res << std::endl;
    }
    
}

void TPZPoroElastoPlasticAnalysis::UpdateState(){
    m_simulation_data->SetTransferCurrentToLastQ(true);
    AcceptPseudoTimeStepSolution();
    m_simulation_data->SetTransferCurrentToLastQ(false);
}


void TPZPoroElastoPlasticAnalysis::PostProcessTimeStep(std::string & file){
    
    int dim = Mesh()->Dimension();
    int div = 0;
    TPZStack< std::string> vecnames;
    m_post_processor->TransferSolution();
    m_post_processor->DefineGraphMesh(dim,m_var_names,vecnames,file);
    m_post_processor->PostProcess(div,dim);
}

void TPZPoroElastoPlasticAnalysis::AcceptPseudoTimeStepSolution(){

    //m_simulation_data->SetInitialStateQ(true);
    SetUpdateMemmory(true);
    AssembleResidual();  //oioioioio
    SetUpdateMemmory(false);
    
//    m_simulation_data->Set_must_accept_solution_Q(false);
//    bool state = m_simulation_data->IsCurrentStateQ();
//    if (state) {
//        m_simulation_data->Set_must_accept_solution_Q(true);
//        AssembleResidual();
//        m_simulation_data->Set_must_accept_solution_Q(false);
//    }else{
//        m_simulation_data->Set_must_accept_solution_Q(true);
//        AssembleResidual();
//        m_simulation_data->Set_must_accept_solution_Q(false);
//    }
}


void TPZPoroElastoPlasticAnalysis::LoadCurrentState(){
    LoadSolution(m_X_n);
    //DebugStop();
}

void TPZPoroElastoPlasticAnalysis::LoadLastState(){
    LoadSolution(m_X);
    DebugStop();
}


void TPZPoroElastoPlasticAnalysis::SetUpdateMemmory(bool accept_solution_Q){

    m_simulation_data->Set_must_accept_solution_Q(accept_solution_Q);
    std::map<int, TPZMaterial * >::iterator mit;
    std::map<int, TPZMaterial *> & refMatVec = Mesh()->MaterialVec();

    TPZMatWithMem<TPZMemoryDFN> * material_with_memory;
    for(mit=refMatVec.begin(); mit!= refMatVec.end(); mit++)
    {
        material_with_memory = dynamic_cast<TPZMatWithMem<TPZMemoryDFN> *>( mit->second );
        if(material_with_memory)
        {
            material_with_memory->SetUpdateMem(accept_solution_Q);
        }
    }
}



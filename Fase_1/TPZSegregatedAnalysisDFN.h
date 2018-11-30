//
//  TPZSegregatedAnalysis.h
//
//  Created by Omar Durán on 9/13/18.
//

#ifndef TPZSegregatedAnalysisDFN_h
#define TPZSegregatedAnalysisDFN_h

#include <stdio.h>
#include "TPZMatWithMem.h"
#include "TPZPoroElastoPlasticAnalysis.h"
#include "TPZPoroElastoPlasticDFN.h"
#include "TPZMonoPhasicMemoryDFN.h"
#include "TPZDarcyAnalysis.h"
#include "TPZDarcy2DMaterialMem.h"
#include "TPZMemoryDFN.h"
#include "TPZMemoryFracDFN.h"
#include "TPZElasticCriterion.h"


class TPZSegregatedAnalysisDFN {
    
private:
    
    /// Pointer to Simulation data object
    TPZSimulationData * m_simulation_data;
    
    /// Pointer to geomechanic analysis object
    TPZPoroElastoPlasticAnalysis * m_elastoplast_analysis;
    
    /// Pointer to reservoir analysis object
    TPZDarcyAnalysis * m_darcy_analysis;
    
public:
    
    /// Default constructor
    TPZSegregatedAnalysisDFN();
    
    /// Destructor
    ~TPZSegregatedAnalysisDFN();
    
    /// Copy constructor
    TPZSegregatedAnalysisDFN(const TPZSegregatedAnalysisDFN & other);
    
    /// Set the pointer of Simulation data object
    void SetSimulationData(TPZSimulationData * simulation_data){
        m_simulation_data = simulation_data;
    }

    /// Attach materials with the same share pointer (For fractures linking)
    void ApplyFracMemoryLink(int frac_matid);
    
    /// Attach materials with the same share pointer (For now just volumeric linking)
    void ApplyMemoryLink(int matid);
    
    /// Configurate internal analysis objects and linked them through the memory shared pointer
    void ConfigurateAnalysis(DecomposeType decompose_geo, DecomposeType decompose_res, TPZSimulationData * simulation_data, TPZCompMesh * cmesh_geomechanics, TPZCompMesh * cmesh_reservoir, TPZManVector<TPZCompMesh * , 2> & mesh_vec, TPZStack<std::string> & post_pro_var_res, TPZStack<std::string> & post_pro_var_geo);
    
    /// Execute the evolution for a single time step
    void ExecuteOneTimeStep();
    
    /// Post-processing the variables for a single time step
    void PostProcessTimeStep(std::string & geo_file, std::string & res_file);
    
    /// Execute the transient evolution using Fixed Stress Split Iteration
    void ExecuteTimeEvolution();

    // Update parameters
    void SetInitialParameters();
    
    // Update parameters
    void UpdateParameters();
    
    /// Update solution state x = x_n
    void UpdateState();
    
    // Adjust integration orders
    void AdjustIntegrationOrder(TPZCompMesh * cmesh_o, TPZCompMesh * cmesh_d);

    // Adjust fracture integration orders
    void AdjustFractureIntegrationOrder(TPZCompMesh * cmesh_o, TPZCompMesh * cmesh_d);
};

#endif /* TPZSegregatedAnalysis.h */

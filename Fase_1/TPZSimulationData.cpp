//
//  TPZSimulationData.cpp
//  PZ
//
//  Created by Omar on 8/28/16.
//
//

#include "TPZSimulationData.h"

TPZSimulationData::TPZSimulationData()
{
    m_h_level           = 0;
    m_elasticity_order  = 1;
    m_darcy_order       = 1;
    m_dimesion          = 0;
    m_n_iterations       = 0;
    m_epsilon_res       = 0;
    m_epsilon_cor       = 0;
    m_n_threads         = 0;
    m_geometry_file     = "";
    m_vtk_file          = "";
    m_vtk_resolution    = 0;
    m_volumetric_material_id.resize(0);
    m_elasticity_ID = 0;
    m_darcy_ID = 0;
    m_is_initial_state_Q  = false;
    m_is_current_state_Q  = false;
    m_must_accept_solution_Q = false;
    m_is_dual_formulation_Q = true;
}

TPZSimulationData::~TPZSimulationData()
{
    
}

TPZSimulationData::TPZSimulationData(const TPZSimulationData & other)
{
    m_h_level                           = other.m_h_level;
    m_elasticity_order                  = other.m_elasticity_order;
    m_darcy_order                       = other.m_darcy_order;
    m_dimesion                          = other.m_dimesion;
    m_n_iterations                      = other.m_n_iterations;
    m_epsilon_res                       = other.m_epsilon_res;
    m_epsilon_cor                       = other.m_epsilon_cor;
    m_n_threads                         = other.m_n_threads;
    m_geometry_file                     = other.m_geometry_file;
    m_vtk_file                          = other.m_vtk_file;
    m_vtk_resolution                    = other.m_vtk_resolution;
    m_volumetric_material_id            = other.m_volumetric_material_id;
    m_elasticity_ID                     = other.m_elasticity_ID;
    m_darcy_ID                          = other.m_darcy_ID;
    m_is_initial_state_Q                = other.m_is_initial_state_Q;
    m_is_current_state_Q                = other.m_is_current_state_Q;
    m_mat_ids                           = other.m_mat_ids;
    m_is_dual_formulation_Q             = other.m_is_dual_formulation_Q;
    m_must_accept_solution_Q            = other.m_must_accept_solution_Q;
    m_transfer_current_to_last_solution_Q            = other.m_transfer_current_to_last_solution_Q;
    m_transfer_current_to_last_solution_Q = false;
}

TPZSimulationData & TPZSimulationData::operator=(const TPZSimulationData &other)
{
    if (this != & other) // prevent self-assignment
    {
        m_h_level                           = other.m_h_level;
        m_elasticity_order                  = other.m_elasticity_order;
        m_darcy_order                       = other.m_darcy_order;
        m_dimesion                          = other.m_dimesion;
        m_n_iterations                      = other.m_n_iterations;
        m_epsilon_res                       = other.m_epsilon_res;
        m_epsilon_cor                       = other.m_epsilon_cor;
        m_n_threads                         = other.m_n_threads;
        m_geometry_file                     = other.m_geometry_file;
        m_vtk_file                          = other.m_vtk_file;
        m_vtk_resolution                    = other.m_vtk_resolution;
        m_volumetric_material_id            = other.m_volumetric_material_id;
        m_elasticity_ID                     = other.m_elasticity_ID;
        m_darcy_ID                          = other.m_darcy_ID;
        m_is_initial_state_Q                = other.m_is_initial_state_Q;
        m_is_current_state_Q                = other.m_is_current_state_Q;
        m_mat_ids                           = other.m_mat_ids;
        m_must_accept_solution_Q            = other.m_must_accept_solution_Q;
        m_transfer_current_to_last_solution_Q            = other.m_transfer_current_to_last_solution_Q;
    }
    return *this;
}

void TPZSimulationData::Print()
{
    
    std::cout << " TPZSimulationData class members : " << std::endl;
    std::cout << std::endl;
    std::cout << " m_h_level = " << m_h_level << std::endl;
    std::cout << " m_elasticity_order = " << m_elasticity_order << std::endl;
    std::cout << " m_elasticity_matid = " << m_elasticity_ID << std::endl;
    std::cout << " m_darcy_order = " << m_darcy_order << std::endl;
    std::cout << " m_darcy_matid = " << m_darcy_ID << std::endl;
    std::cout << " m_dimesion = " << m_dimesion << std::endl;
    std::cout << " m_n_iterations = " << m_n_iterations << std::endl;
    std::cout << " m_epsilon_res = " << m_epsilon_res << std::endl;
    std::cout << " m_epsilon_cor = " << m_epsilon_cor << std::endl;
    std::cout << " m_n_threads = " << m_n_threads << std::endl;
    std::cout << " m_geometry_file = " << m_geometry_file << std::endl;
    std::cout << " m_vtk_file = " << m_vtk_file << std::endl;
    std::cout << " m_vtk_resolution = " << m_vtk_resolution << std::endl;
    std::cout << " m_volumetric_material_id = " << &m_volumetric_material_id << std::endl;
    std::cout << std::endl;
    
}


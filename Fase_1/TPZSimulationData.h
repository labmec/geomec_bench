//
//  TPZSimulationData.h
//  PZ
//
//  Created by Omar on 8/28/18.
//
//

#ifndef TPZSimulationData_h
#define TPZSimulationData_h

#include <stdio.h>
#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzstack.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include "TPZGmshReader.h"
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzcheckgeom.h"


/** @brief Object conatining several kind of informations being used anytime and anywhere */
class TPZSimulationData
{
    
protected:
    
    /** @brief Spatial refinemenet level */
    int m_h_level;
    
    /** @brief Polynomial order for elasticity component */
    int m_elasticity_order;
    
    /** @brief Material ID related to Elasticity material */
    int m_elasticity_ID;
    
    /** @brief Polynomial order for diffusion component */
    int m_darcy_order;
    
    /** @brief Material ID related to Darcy material */
    int m_darcy_ID;
    
    /** @brief Physical dimension of the domain */
    int m_dimesion;
    
    /** @brief Number of iteration */
    int m_n_iterations;
    
    /** @brief Residue overal tolerance */
    REAL m_epsilon_res;
    
    /** @brief Correction overal tolerance */
    REAL m_epsilon_cor;
    
    /** @brief Number of thread */
    int m_n_threads;

    bool m_is_initial_stress_Q;
    
    bool m_is_initial_state_Q;
    
    bool m_is_current_state_Q;
    
    /** @brief Verify if the problem is monophasic */
    bool m_is_mono_Q;
    
    /** @brief Directive that states if the last memory solution is being transferred to the current memory solution */
    bool m_transfer_current_to_last_solution_Q;
    
    /** @brief Directive that states if the current solution must be accepted inside the memory  */
    bool m_must_accept_solution_Q;
    
    /** @brief Directive that states the use of dual (true) or pirmal (false) formulation for monophacic flow  */
    bool m_is_dual_formulation_Q;
    
    /** @brief Verify if fracture has to inserted */
    bool m_insert_fractures_Q;
    
    /** @brief Name for the Gmsh geometry file being used */
    std::string m_geometry_file;
    
    /** @brief Name for the vtk files being postprocessed */
    std::string m_vtk_file;
    
    /** @brief Number of vtk resolution during postprocessing */
    int m_vtk_resolution;
    
    /** @brief Vector that storages only volumetric material identifiers (higher dimension elements) */
    std::vector<int> m_volumetric_material_id;
    
    /** @brief Vector that storages only fracture material identifiers (higher dimension elements) */
    std::vector<int> m_fracture_material_id;
    
    /** @brief Vector taht storages interface identifiers */
    std::vector<int> m_interface_id;

    /** @brief Vector taht storages left interface identifiers */
    std::vector<int> m_interfaceLeft_id;

    /** @brief Vector taht storages right interface identifiers */
    std::vector<int> m_interfaceRight_id;
    
    /** @brief Material and boundaries identifiers sorted per region */
    TPZManVector<std::pair<int, TPZManVector<int,12>>,12> m_mat_ids;

    /** @brief Young modulus */
    REAL m_Eyoung;

    /** @brief Poisson coeficient */
    REAL m_poisson;

    /** @brief Biot coeficient */
    REAL m_alpha;
    
    /** @brief Initial stress state. Firts simulation scenario*/
    TPZFNMatrix<9,REAL> m_Stress0;

    /** @brief Initial volumetric stress state, for each integration point. Firts simulation scenario*/
    std::vector<REAL> m_Stress_Vol0;
    
    /** @brief Initial porous media permeability tensor. Coeficient which multiplies the gradient operator*/
    TPZFNMatrix<9,REAL> m_TensorK0;
    
    /** @brief Initial permeability. Coeficient which multiplies the gradient operator*/
    std::map<REAL, REAL> m_k0;
    
    /** @brief Fracture orientation. Define fracture direction*/
    std::map<REAL, REAL> m_fracOrient;
    
    /** @brief Initial porosity.*/
    REAL m_phi0;
    
    /** Fracture max closure for each fracture mat id */
    std::map<REAL,REAL> m_Vm;

    /** Fracture initial opening for each fracture mat id */
    std::map<REAL,REAL> m_a0;
    
    /** Fracture normal stiffiness for each fracture mat id */
    REAL m_Kni;
    
    
    
public:
    
    /** @brief default constructor */
    TPZSimulationData();
    
    /** @brief default constructor */
    TPZSimulationData(const TPZSimulationData & other);
    
    /** @brief default constructor */
    TPZSimulationData &operator=(const TPZSimulationData &other);
    
    /** @brief destructor */
    ~TPZSimulationData();
    
    /** @brief Print object attributes */
    void Print();
    
    /** @brief Set the directive that states if the current memory solution is being transferred to the last memory solution */
    void SetTransferCurrentToLastQ(bool transfer_current_to_last_solution_Q) { m_transfer_current_to_last_solution_Q = transfer_current_to_last_solution_Q; }
    
    /** @brief Get the directive that states if the current memory solution is being transferred to the last memory solution */
    bool GetTransferCurrentToLastQ() { return m_transfer_current_to_last_solution_Q; }
    
    /// Access methods
    
    /** @brief Get the material and boundaries identifiers sorted per region */
    TPZManVector<std::pair<int, TPZManVector<int,12>>,12> & MaterialIds() { return m_mat_ids; }

    /** @brief Set initial state */
    void SetInitialStressQ(bool stressstate) {
        m_is_initial_stress_Q = stressstate;
    }
    
    /** @brief Get initial state */
    bool IsInitialStressQ() {
        return m_is_initial_stress_Q;
    }
    
    
    /** @brief Set initial state */
    void SetInitialStateQ(bool state) {
        m_is_initial_state_Q = state;
    }
    
    /** @brief Get initial state */
    bool IsInitialStateQ() {
        return m_is_initial_state_Q;
    }
    
    /** @brief Set Monophasic Problem */
    void SetMonoPhasicQ(bool is_mono_Q) {
        m_is_mono_Q = is_mono_Q;
    }
    
    /** @brief Verify if is a monophasic problem */
    bool IsMonoPhasicQ() {
        return m_is_mono_Q;
    }
    
    
    /** @brief Set current time state */
    void SetCurrentStateQ(bool state) { m_is_current_state_Q = state; }
    
    /** @brief Get current time state */
    bool IsCurrentStateQ() {return m_is_current_state_Q;}
    
    /** @brief Set the directive that states if the current solution must be accepted inside the memory  */
    void Set_must_accept_solution_Q(bool must_accept_solution_Q){
        m_must_accept_solution_Q = must_accept_solution_Q;
    }
    
    /** @brief Get the directive that states if the current solution must be accepted inside the memory  */
    bool Get_must_accept_solution_Q() { return m_must_accept_solution_Q; }
    
    /** @brief Get the the use of dual (true) or pirmal (false) formulation for monophacic flow  */
    bool Get_is_dual_formulation_Q() { return m_is_dual_formulation_Q; }
    
    /** @brief Verify if fracture has to be inserted */
    void Set_insert_fractures_Q(bool insert_fractures_Q){
        m_insert_fractures_Q = insert_fractures_Q;
    }
    
    /** @brief Verify if fracture has to be inserted   */
    bool Get_insert_fractures_Q() { return m_insert_fractures_Q; }
    
    /** @brief Set the spatial refinemenet level */
    void Set_h_level(int h_level){
        m_h_level = h_level;
    }
    
    /** @brief Get the spatial refinemenet level */
    int Get_h_level(){
        return m_h_level;
    }
    
    /** @brief Set the polynomial order for the elasticity approximation */
    void Set_elasticity_order(int elasticity_order){
        m_elasticity_order = elasticity_order;
    }
    /** @brief Set the material ID related to Elasticity material */
    void Set_elasticity_matid(int elasticity_matid){
        m_elasticity_ID = elasticity_matid;
    }
    
    /** @brief Get the material ID related to Elasticity material */
    int Get_elasticity_matid(){
        return m_elasticity_ID;
    }
    
    /** @brief Get the polynomial order for the elasticity approximation  */
    int Get_elasticity_order(){
        return m_elasticity_order;
    }
    
    /** @brief Set the polynomial order for the Darcy's approximation */
    void Set_darcy_order(int darcy_order){
        m_darcy_order = darcy_order;
    }
    
    /** @brief Set the material ID related to Darcy's material */
    void Set_darcy_matid(int darcy_matid){
        m_darcy_ID = darcy_matid;
    }
    
    /** @brief Get the material ID related to Darcy's material */
    int Get_darcy_matid(){
        return m_darcy_ID;
    }
    
    /** @brief Get the polynomial order for the Darcy's approximation  */
    int Get_darcy_order(){
        return m_darcy_order;
    }
    
    /** @brief Set the problem dimension */
    void Set_dimesion(int dimesion){
        m_dimesion = dimesion;
    }
    
    /** @brief Get the problem dimension  */
    int Get_dimesion(){
        return m_dimesion;
    }
    
    /** @brief Set Newton iterations */
    void Set_n_iterations(int n_iterations){
        m_n_iterations = n_iterations;
    }
    
    /** @brief Get Newton iterations  */
    REAL Get_n_iterations(){
        return m_n_iterations;
    }
    
    /** @brief Set residue tolerance */
    void Set_epsilon_res(REAL epsilon_res){
        m_epsilon_res = epsilon_res;
    }
    
    /** @brief Get residue tolerance  */
    REAL Get_epsilon_res(){
        return m_epsilon_res;
    }
    
    /** @brief Set correction tolerance */
    void Set_epsilon_cor(REAL epsilon_cor){
        m_epsilon_cor = epsilon_cor;
    }
    
    /** @brief Get correction tolerance  */
    REAL Get_epsilon_cor(){
        return m_epsilon_cor;
    }
    
    /** @brief Get number of threads being used  */
    int Get_n_threads(){
        return m_n_threads;
    }
    
    /** @brief Set number of threads being used */
    void Set_n_threads(int n_threads){
        m_n_threads = n_threads;
    }
    
    /** @brief Get name for the Gmsh geometry file being used  */
    std::string Get_geometry_file(){
        return m_geometry_file;
    }
    
    /** @brief Set name for the Gmsh geometry file being used */
    void Set_n_threads(std::string geometry_file){
        m_geometry_file = geometry_file;
    }
    
    /** @brief Get name for the vtk files being postprocessed  */
    std::string Get_vtk_file(){
        return m_vtk_file;
    }
    
    /** @brief Set name for the vtk files being postprocessed */
    void Set_vtk_file(std::string vtk_file){
        m_vtk_file = vtk_file;
    }
    
    /** @brief Get number of vtk resolution during postprocessing  */
    int Get_vtk_resolution(){
        return m_vtk_resolution;
    }
    
    /** @brief Set number of vtk resolution during postprocessing */
    void Set_vtk_resolution(int vtk_resolution){
        m_vtk_resolution = vtk_resolution;
    }
    
    /** @brief Get the vector that storages only volumetric material identifiers (higher dimension elements)  */
    std::vector<int> & Get_volumetric_material_id(){
        return m_volumetric_material_id;
    }
    
    /** @brief Set the vector that storages only volumetric material identifiers (higher dimension elements) */
    void Set_volumetric_material_id(std::vector<int> & volumetric_material_id){
        m_volumetric_material_id = volumetric_material_id;
    }

    /** @brief Get the vector that storages only fracture material identifiers (higher dimension elements)  */
    std::vector<int> & Get_fracture_material_id(){
        return m_fracture_material_id;
    }
    
    /** @brief Set the vector that storages only fracture material identifiers (higher dimension elements) */
    void Set_fracture_material_id(std::vector<int> & fracture_material_id){
        m_fracture_material_id = fracture_material_id;
    }
    
    /** @brief Get the vector that storages interface identifiers  */
    std::vector<int> & Get_interface_id(){
        return m_interface_id;
    }
    
    /** @brief Set the vector that storages interface identifiers  */
    void Set_interface_id(std::vector<int> & interface_id){
        m_interface_id = interface_id;
    }

    /** @brief Get the vector that storages left interface identifiers  */
    std::vector<int> & Get_interfaceLeft_id(){
        return m_interfaceLeft_id;
    }
    
    /** @brief Set the vector that storages left interface identifiers  */
    void Set_interfaceLeft_id(std::vector<int> & interfaceLeft_id){
        m_interfaceLeft_id = interfaceLeft_id;
    }
    
    /** @brief Get the vector that storages interface identifiers  */
    std::vector<int> & Get_interfaceRight_id(){
        return m_interfaceRight_id;
    }
    
    /** @brief Set the vector that storages interface identifiers  */
    void Set_interfaceRight_id(std::vector<int> & interfaceRight_id){
        m_interfaceRight_id = interfaceRight_id;
    }
    /** @brief Get Young modulus of material  */
    REAL & Get_Eyoung(){
        return m_Eyoung;
    }
    
    /** @brief Set Young modulus of material  */
    void Set_Eyoung(REAL Eyoung){
        m_Eyoung = Eyoung;
    }
    
    /** @brief Get Poisson coeficient of material  */
    REAL & Get_Poisson(){
        return m_poisson;
    }
    
    /** @brief Set Poisson coeficient of material  */
    void Set_Poisson(REAL poisson){
        m_poisson = poisson;
    }
    
    /** @brief Get Biot coeficient of material  */
    REAL & Get_Biot(){
        return m_alpha;
    }
    
    /** @brief Set Biot coeficient of material  */
    void Set_Biot(REAL alpha){
        m_alpha = alpha;
    }
    

    /** @brief Set the initial volumetric stress per integration point index */
    void Set_Stress_Vol0(std::vector<REAL> Stress0){
        m_Stress_Vol0 = Stress0;
    }
    
    /** @brief Get the initial volumetric stress per integration point index */
    std::vector<REAL> Get_Stress_Vol0(){
        return m_Stress_Vol0;
    }
    
    /** @brief Set the initial stress state tensor */
    void Set_Stress0(TPZFNMatrix<9,REAL> Stress0){
        m_Stress0 = Stress0;
    }
    
    /** @brief Get the initial stress state tensor */
    TPZFNMatrix<9,REAL> Get_Stress0(){
        return m_Stress0;
    }
    
    /** @brief Set the initial Porous media permeability tensor */
    void Set_PermeabilityTensor_0(TPZFNMatrix<9,REAL> TensorK0){
        m_TensorK0 = TensorK0;
    }
    
    /** @brief Get the initial Porous media permeability tensor */
    TPZFNMatrix<9,REAL> Get_PermeabilityTensor_0(){
        return m_TensorK0;
    }

    /** @brief Set the initial permeability */
    void Set_Permeability_0(std::map<REAL, REAL> k0){
        m_k0 = k0;
    }
    
    /** @brief Get the initial permeability */
    std::map<REAL, REAL> Get_Permeability_0(){
        return m_k0;
    }

    /** @brief Set fracture orientation */
    void Set_FractureOrient(std::map<REAL, REAL> orient){
        m_fracOrient = orient;
    }
    
    /** @brief Get fracture orientation */
    std::map<REAL, REAL> Get_FractureOrient(){
        return m_fracOrient;
    }
    
    /// Set lagrangian porosity at intial REAL
    void Set_Porosity_0(REAL phi_0)
    {
        m_phi0 = phi_0;
    }
    
    /// Get lagrangian porosity at intial REAL
    REAL Get_Porosity_0()
    {
        return m_phi0;
    }
    
    // Set max closure for each fracture
    void Set_Vm(std::map<REAL, REAL> Vm_frac){
        m_Vm = Vm_frac;
    }
    
    // Get max closure for each fracture
    std::map<REAL,REAL> Get_Vm(){
        return m_Vm;
    }
    
    // Set initial opening for each fracture
    void Set_a0(std::map<REAL, REAL> a0_frac){
        m_a0 = a0_frac;
    }

    // Get initial opening for each fracture
    std::map<REAL,REAL> Get_a0(){
        return m_a0;
    }

    // Set normal stiffness for each fracture
    void Set_Kni(REAL kni_frac){
        m_Kni = kni_frac;
    }

    // Get normal stiffness for each fracture
    REAL Get_Kni(){
        return m_Kni;
    }
    
    
};

#endif /* TPZSimulationData_h */

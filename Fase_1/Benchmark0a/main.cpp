
#include "pzgmesh.h"
#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>

#include "HidraulicoMonofasico2D.h"
#include "HidraulicoMonofasicoElastico.h"
//#include "FraturaElastico.h"

#include "pzlog.h"
//------------------Benchmarks for Geomec------------------------

using namespace std;

/// Main dedicated to execute each scenario for benchmark_1
int main(int argc, char *argv[])
{
//    std::string log_cfg_file = "/Users/pablocarvalho/Documents/GitHub/geomec_bench/Fase_1/benchmark.cfg";
//    InitializePZLOG(log_cfg_file);
//    HidraulicoMonofasicoElastico scenario0a;
//    //StopError();
//    scenario0a.Run(2);
    
    InitializePZLOG();
    HidraulicoMonofasicoElastico scenario0a;
    //StopError();
    
    // Caso inicial DeltaP = 0
    TPZSimulationData *simulation_data = new TPZSimulationData;
    simulation_data->SetInitialStressQ(true);
    //simulation_data->SetMonoPhasicQ(true);
    scenario0a.SetSimulationData(simulation_data);
    scenario0a.Run(2);
    
    // Demais casos
    int N_cases = 9;
    TPZVec<TPZSimulationData *> sim_data_vec(N_cases);
    TPZVec<HidraulicoMonofasicoElastico *> scenarios_vec(N_cases);
    REAL DeltaP = 5;
    
    for (int i_case = 0; i_case < N_cases; i_case++) {
        sim_data_vec[i_case] = new TPZSimulationData;
        sim_data_vec[i_case]->Set_Stress_Vol0(simulation_data->Get_Stress_Vol0());
        sim_data_vec[i_case]->SetInitialStressQ(false);
        scenarios_vec[i_case] = new HidraulicoMonofasicoElastico;
        scenarios_vec[i_case]->SetDeltaP(DeltaP);
        scenarios_vec[i_case]->SetSimulationData(sim_data_vec[i_case]);
        scenarios_vec[i_case]->Run(2);
        DeltaP = DeltaP+5.;
    }
    
    
    
//    MonofasicoElastico cenario0a;
//    cenario0a.Run(2);
    return 0;
}

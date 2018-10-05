
#include "pzgmesh.h"
#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>

#include "HidraulicoMonofasico2D.h"
#include "HidraulicoMonofasicoElastico.h"
#include "FraturaElastico.h"

#include "pzlog.h"
//------------------Benchmarks for Geomec------------------------

using namespace std;

/// Main dedicated to execute each scenario for benchmark_1
int main(int argc, char *argv[])
{
    std::string log_cfg_file = "/Users/pablocarvalho/Documents/GitHub/geomec_bench/Fase_1/benchmark.cfg";
    InitializePZLOG(log_cfg_file);
    HidraulicoMonofasicoElastico scenario0a;
    //StopError();
    scenario0a.Run(2);
    
    
//    MonofasicoElastico cenario0a;
//    cenario0a.Run(2);
    return 0;
}

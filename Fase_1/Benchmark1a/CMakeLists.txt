

file(GLOB commonHeaders ../*.h)
file(GLOB commonSources ../*.cpp)

file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(Mymain main.cpp ${sources} ${headers} ${commonHeaders} ${commonSources})

add_executable(Benchmark1a ${Mymain})

target_compile_definitions(Benchmark1a PUBLIC BENCHMARK_SOURCE_DIR=${PROJECT_SOURCE_DIR})

target_link_libraries(Benchmark1a pz )

configure_file(gmsh/GeometryBench.msh GeometryBench.msh COPYONLY)
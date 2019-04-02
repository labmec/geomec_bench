########### next target ###############

add_executable(FracFlow FracFlow.cpp TPZMatLaplacianHybrid.h TPZMatLaplacianHybrid.cpp meshgen.h meshgen.cpp)

target_link_libraries(FracFlow pz ${CGAL})

add_executable(FracPre FracPre.cpp TPZFracSet.h TPZFracSet.cpp)

target_link_libraries(FracPre pz ${CGAL})

add_executable(FracFlowGmsh FracFlowGmsh.cpp TPZMatLaplacianHybrid.h TPZMatLaplacianHybrid.cpp TPZFracSet.cpp TPZFracSet.h TPZFracSimulation.h TPZFracSimulation.cpp)

target_link_libraries(FracFlowGmsh pz ${CGAL})

add_executable(FractureIntersection TPZFractureIntersection.cpp TPZPointCloud.h TPZPointCloud.cpp FractureIntersectionConfig.h FractureIntersectionConfig.cpp TPZElementIntersect.h TPZElementIntersect.cpp
    VTKIntersect.h VTKIntersect.cpp)

target_link_libraries(FractureIntersection pz ${CGAL} )

file(COPY BoxReservoir.geo DESTINATION .)
file(COPY cube.geo DESTINATION .)
file(COPY plane.geo DESTINATION .)

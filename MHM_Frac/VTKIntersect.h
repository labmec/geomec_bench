//
//  VTKIntersect.hpp
//  Benchmark0a
//
//  Created by Philippe Devloo on 21/09/18.
//

#ifndef VTKIntersect_hpp
#define VTKIntersect_hpp

#include <stdio.h>

#ifdef USING_VTK

#include <vtkActor.h>;

void GenerateVTKInput();

void VTKWindow();

void GetBooleanOperationActor( double x, int operation );

vtkActor * GetStructuredGrid();

#endif

#endif /* VTKIntersect_hpp */

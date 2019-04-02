//
//  VTKIntersect.cpp
//  Benchmark0a
//
//  Created by Philippe Devloo on 21/09/18.
//

#include "VTKIntersect.h"

/*=========================================================================
 
 Program:   Visualization Toolkit
 Module:    TestBooleanOperationPolyDataFilter.cxx
 
 Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 All rights reserved.
 See Copyright.txt or http://www.kitware.com/Copyright.htm for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.
 =========================================================================*/
#ifdef USING_VTK

#include <vtkActor.h>
#include <vtkAppendPolyData.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDistancePolyDataFilter.h>
#include <vtkIntersectionPolyDataFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataReader.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkReverseSense.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkPlaneSource.h>
#include <vtkCubeSource.h>
#include <vtkConeSource.h>
#include <vtkLineSource.h>

#include <vtkThreshold.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkGeometryFilter.h>
#include <vtkDataSetMapper.h>
#include "vtkCamera.h"
#include <vtkStructuredGrid.h>

#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"

vtkActor * GetStructuredGrid()
{
    vtkSmartPointer<vtkStructuredGrid> structuredGrid =
    vtkSmartPointer<vtkStructuredGrid>::New();
    
    vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();
    double x, y, z;
    
    x = 0.0;
    y = 0.0;
    z = 0.0;
    
    for(unsigned int k = 0; k < 2; k++)
    {
        z += 2.0;
        for(unsigned int j = 0; j < 3; j++)
        {
            y += 1.0;
            for(unsigned int i = 0; i < 2; i++)
            {
                x += .5;
                points->InsertNextPoint(x, y, z);
            }
        }
    }
    
    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(2,3,2);
    structuredGrid->SetPoints(points);
    
    int* dims = structuredGrid->GetDimensions();
    
    // Retrieve the entries from the grid and print them to the screen
    unsigned int counter = 0;
    
    for (int k = 0; k < dims[2]; k++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            for (int i = 0; i < dims[0]; i++)
            {
                double p[3];
                structuredGrid->GetPoint(counter, p);
                
                double pNew[3];
                structuredGrid->GetPoint(i, j, k, pNew);
                
                std::cout << "P   : "
                << p[0] << " "
                << p[1] << " "
                << p[2] << std::endl;
                std::cout << "PNew: "
                << pNew[0] << " "
                << pNew[1] << " "
                << pNew[2] << std::endl;
                
                counter++;
            }
        }
    }
    
    
   
    // Create a mapper and actor
    vtkSmartPointer<vtkDataSetMapper> mapper =
    vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputData(structuredGrid);
    
    
    vtkActor * actor = vtkActor::New();
    actor->SetMapper(mapper);
    
    return actor;
}

void GenerateVTKInput()
{
    {
        TPZGmshReader gmsh;
        gmsh.GetDimNamePhysical()[3]["domain"] = 1;

        TPZGeoMesh *gmesh = 0;
#ifdef MACOSX
        gmesh = gmsh.GeometricGmshMesh("../cube.msh");
#else
        gmesh = gmsh.GeometricGmshMesh("cube.msh");
#endif
        std::ofstream out("cube.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        delete gmesh;
    }
    {
        TPZGmshReader gmsh;
        gmsh.GetDimNamePhysical()[2]["domain"] = 1;

        TPZGeoMesh *gmesh = 0;
#ifdef MACOSX
        gmesh = gmsh.GeometricGmshMesh("../plane.msh");
#else
        gmesh = gmsh.GeometricGmshMesh("plane.msh");
#endif
        std::ofstream out("plane.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        delete gmesh;
    }

}

void GetBooleanOperationActor( double x, int operation )
{
    operation = vtkBooleanOperationPolyDataFilter::VTK_UNION;

    double centerSeparation = 0.15;
    vtkSmartPointer<vtkSphereSource> sphere1 =
    vtkSmartPointer<vtkSphereSource>::New();
    //    sphere1->SetCenter(-centerSeparation + x, 0.0, 0.0);
    sphere1->SetCenter(1.0, 1.0, 0.0);
    
    vtkSmartPointer<vtkSphereSource> sphere2 =
    vtkSmartPointer<vtkSphereSource>::New();
    sphere2->SetCenter(  1.5, 1.0, 0.0);
    
    // Create a plane
    vtkSmartPointer<vtkPlaneSource> planeSource =
    vtkSmartPointer<vtkPlaneSource>::New();
    planeSource->SetCenter(0.0, 0.0, 0.0);
    planeSource->SetNormal(0.0, 0.0, 1.0);
    planeSource->SetResolution(12, 30);
    planeSource->Update();
    
    vtkSmartPointer<vtkPlaneSource> planeSource2 =
    vtkSmartPointer<vtkPlaneSource>::New();
    planeSource2->SetCenter(0.0, 0.0, 0.0);
    planeSource2->SetNormal(1.0, 0.0, 1.0);
    planeSource2->SetResolution(10, 10);
    planeSource2->Update();

    vtkSmartPointer<vtkLineSource> lineSource =
    vtkSmartPointer<vtkLineSource>::New();
    
    lineSource->SetPoint1(-2, 0, 0);
    lineSource->SetPoint2(2, 0, 0);
    lineSource->Update();

    // Create a plane
    vtkSmartPointer<vtkCubeSource> cubeSource =
    vtkSmartPointer<vtkCubeSource>::New();
    cubeSource->SetCenter(0.0, 0.0, 0.0);
    cubeSource->SetXLength(2.0);
    cubeSource->SetYLength(2.0);
    cubeSource->SetZLength(2.0);
    
    cubeSource->Update();
//    vtkPolyData* plane = planeSource->GetOutput();

    
    vtkSmartPointer<vtkIntersectionPolyDataFilter> intersection =
    vtkSmartPointer<vtkIntersectionPolyDataFilter>::New();
//    intersection->SetInputData(0, planeSource2->GetOutput());
//    intersection->SetInputData(1, planeSource->GetOutput());
//    intersection->CheckInputOn();
//    intersection->Update();
//    vtkIndent ind(0);
//    intersection->PrintSelf(std::cout, ind);
    
    /*
    vtkSmartPointer<vtkDistancePolyDataFilter> distance =
    vtkSmartPointer<vtkDistancePolyDataFilter>::New();
    distance->SetInputConnection( 0, intersection->GetOutputPort( 1 ) );
    distance->SetInputConnection( 1, intersection->GetOutputPort( 2 ) );
    
    vtkSmartPointer<vtkThreshold> thresh1 =
    vtkSmartPointer<vtkThreshold>::New();
    thresh1->AllScalarsOn();
    thresh1->SetInputArrayToProcess
    ( 0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "Distance" );
    thresh1->SetInputConnection( distance->GetOutputPort( 0 ) );
    
    vtkSmartPointer<vtkThreshold> thresh2 =
    vtkSmartPointer<vtkThreshold>::New();
    thresh2->AllScalarsOn();
    thresh2->SetInputArrayToProcess
    ( 0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "Distance" );
    thresh2->SetInputConnection( distance->GetOutputPort( 1 ) );
    
    if ( operation == vtkBooleanOperationPolyDataFilter::VTK_UNION )
    {
        thresh1->ThresholdByUpper( 0.0 );
        thresh2->ThresholdByUpper( 0.0 );
    }
    else if ( operation == vtkBooleanOperationPolyDataFilter::VTK_INTERSECTION )
    {
        thresh1->ThresholdByLower( 0.0 );
        thresh2->ThresholdByLower( 0.0 );
    }
    else // Difference
    {
        thresh1->ThresholdByUpper( 0.0 );
        thresh2->ThresholdByLower( 0.0 );
    }
    
    vtkSmartPointer<vtkDataSetSurfaceFilter> surface1 =
    vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surface1->SetInputConnection( thresh1->GetOutputPort() );
    
    vtkSmartPointer<vtkDataSetSurfaceFilter> surface2 =
    vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surface2->SetInputConnection( thresh2->GetOutputPort() );
    
    vtkSmartPointer<vtkReverseSense> reverseSense =
    vtkSmartPointer<vtkReverseSense>::New();
    reverseSense->SetInputConnection( surface2->GetOutputPort() );
    if ( operation == 2 ) // difference
    {
        reverseSense->ReverseCellsOn();
        reverseSense->ReverseNormalsOn();
    }
    */
    vtkSmartPointer<vtkAppendPolyData> appender =
    vtkSmartPointer<vtkAppendPolyData>::New();
//    appender->SetInputConnection( surface1->GetOutputPort() );
//    if ( operation == 2)
//    {
//        appender->AddInputConnection( reverseSense->GetOutputPort() );
//    }
//    else
//    {
//        appender->AddInputConnection( surface2->GetOutputPort() );
//    }
    
//    appender->AddInputData(cubeSource->GetOutput());
//    appender->AddInputData(sphere1->GetOutput());
    vtkSmartPointer<vtkPolyDataMapper> mapperCube =
    vtkSmartPointer<vtkPolyDataMapper>::New();
    mapperCube->SetInputConnection( cubeSource->GetOutputPort() );
    mapperCube->ScalarVisibilityOff();
    vtkActor *actorCube = vtkActor::New();
    actorCube->SetMapper( mapperCube );

    vtkSmartPointer<vtkPolyDataMapper> mapperPlane =
    vtkSmartPointer<vtkPolyDataMapper>::New();
    mapperPlane->SetInputConnection( planeSource->GetOutputPort() );
    mapperPlane->ScalarVisibilityOff();
    vtkActor *actorPlane = vtkActor::New();
    actorPlane->SetMapper( mapperPlane );
    
    vtkSmartPointer<vtkPolyDataMapper> mapperPlane2 =
    vtkSmartPointer<vtkPolyDataMapper>::New();
    mapperPlane2->SetInputConnection( planeSource2->GetOutputPort() );
    mapperPlane2->ScalarVisibilityOff();
    vtkActor *actorPlane2 = vtkActor::New();
    actorPlane2->SetMapper( mapperPlane2 );
  
    vtkSmartPointer<vtkPolyDataMapper> mapperSphere1 =
    vtkSmartPointer<vtkPolyDataMapper>::New();
    mapperSphere1->SetInputConnection( sphere1->GetOutputPort() );
    mapperSphere1->ScalarVisibilityOff();

    
    vtkActor *actorSphere = vtkActor::New();
    actorSphere->SetMapper( mapperSphere1 );
    vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
    //renderer->AddActor(actorCube);
    //renderer->AddActor(actorSphere);
    renderer->AddActor(actorPlane);
    renderer->AddActor(actorPlane2);
    vtkSmartPointer<vtkActor> gridActor;// = vtkSmartPointer<vtkActor>::New();
    gridActor = GetStructuredGrid();
    renderer->AddActor(gridActor);

    vtkSmartPointer<vtkActor> gridActor2;// = vtkSmartPointer<vtkActor>::New();
    gridActor2 = GetStructuredGrid();
    gridActor2->RotateX(4);
    gridActor2->RotateY(8);
    gridActor2->RotateZ(3);
    
    
    vtkDataSet* dataset = gridActor->GetMapper()->GetInput();
    vtkPolyData* pData = vtkPolyData::SafeDownCast(dataset);
    vtkDataSet* dataset2 = gridActor2->GetMapper()->GetInput();
    vtkPolyData* pData2 = vtkPolyData::SafeDownCast(dataset2);

    intersection->SetInputData(0, pData);
    intersection->SetInputData(1, pData2);
    intersection->CheckInputOn();
    intersection->Update();
    vtkIndent ind(0);
    intersection->PrintSelf(std::cout, ind);

    
    
    renderer->AddActor(gridActor2);

    renderer->SetBackground(0.1, 0.2, 0.4);
    // Zoom in a little by accessing the camera and invoking its "Zoom" method.
    renderer->ResetCamera();
    renderer->GetActiveCamera()->Zoom(1.5);
    
    // The render window is the actual GUI window
    // that appears on the computer screen
    vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->SetSize(2000, 2000);
    renderWindow->AddRenderer(renderer);
    
    // The render window interactor captures mouse events
    // and will perform appropriate camera or actor manipulation
    // depending on the nature of the events.
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);
    
    // This starts the event loop and as a side effect causes an initial render.
    renderWindowInteractor->Start();

}

/*
int TestBooleanOperationPolyDataFilter2(int, char *[])
{
    vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
    
    vtkSmartPointer<vtkRenderWindow> renWin =
    vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer( renderer );
    
    vtkSmartPointer<vtkRenderWindowInteractor> renWinInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renWinInteractor->SetRenderWindow( renWin );
    
    vtkActor *unionActor =
    GetBooleanOperationActor( -2.0, vtkBooleanOperationPolyDataFilter::VTK_UNION );
    renderer->AddActor( unionActor );
    unionActor->Delete();
    
    vtkActor *intersectionActor =
    GetBooleanOperationActor(  0.0, vtkBooleanOperationPolyDataFilter::VTK_INTERSECTION );
    renderer->AddActor( intersectionActor );
    intersectionActor->Delete();
    
    vtkActor *differenceActor =
    GetBooleanOperationActor(  2.0, vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE );
    renderer->AddActor( differenceActor );
    differenceActor->Delete();
    
    renWin->Render();
    renWinInteractor->Start();
    
    return EXIT_SUCCESS;
}
*/

#endif

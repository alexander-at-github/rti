#pragma once

#include <boost/core/demangle.hpp>

#include <vtkAbstractArray.h>
#include <vtkCellIterator.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDataSetAttributes.h>
#include <vtkGenericCell.h>
#include <vtkPointData.h>
//#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtksys/SystemTools.hxx>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>

#include "rti/io/i_triangle_reader.hpp"
#include "rti/util/utils.hpp"

// There is an example at https://vtk.org/doc/nightly/html/classvtkCellIterator.html
//
// Good example at https://vtk.org/gitweb?p=VTK.git;a=blob;f=Examples/IO/Cxx/DumpXMLFile.cxx


namespace rti { namespace io { namespace christoph {
  // The parameter Ty is intended to be instantiated as a numeric type.
  template<typename Ty>
  class vtu_triangle_reader : public rti::io::i_triangle_reader<Ty> {
  public:
    vtu_triangle_reader(const std::string& pFilename) :
      rti::io::i_triangle_reader<Ty>(pFilename) {
      std::string extension = vtksys::SystemTools::GetFilenameLastExtension(pFilename);
      auto extensionCondition = std::string {".vtu"};
      if (extension != extensionCondition) {
        std::cerr
          << "Warning: " << boost::core::demangle(typeid(this).name())
          << " may not be able to read the file " << pFilename
          << " because it does not have the " << extensionCondition
          << " extension." << std::endl;
      }
      auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
      reader->SetFileName(pFilename.c_str());
      reader->Update();
      //vtkSmartPointer<vtkUnstructuredGrid> unstructuredgrid = reader->GetOutput(); // type: vtkUnstructuredGrid
      auto unstructuredgrid = reader->GetOutput(); // type: vtkUnstructuredGrid
      auto numPnts = unstructuredgrid->GetNumberOfPoints(); // type: vtkIdType
      auto pointdata = unstructuredgrid->GetPointData();
      // if (pointdata != nullptr) {
      //   std::cout << " contains point data with "
      //             << pointdata->GetNumberOfArrays()
      //             << " arrays." << std::endl;
      //   for (int i = 0; i < pointdata->GetNumberOfArrays(); i++) {
      //     std::cout << "\tArray " << i
      //               << " is named "
      //               << (pointdata->GetArrayName(i) ? pointdata->GetArrayName(i) : "NULL")
      //               << std::endl;
      //   }
      // }
      if (pointdata == nullptr) {
        std::cerr
          << "Warning: "  << boost::core::demangle(typeid(this).name())
          << " could not find data in the file " << pFilename << std::endl;
      }
      // Write points from vtk to this data structure
      this->mPoints.reserve(numPnts);
      for (vtkIdType idx = 0; idx < numPnts; ++idx) {
        double xyz[3]; // 3 dimensions
        unstructuredgrid->GetPoint(idx, xyz);
        this->mPoints.push_back({(Ty) xyz[0], (Ty) xyz[1], (Ty) xyz[2]});
      }

      // Write triangles from the VTK data structure
      auto cellarray = vtkSmartPointer<vtkCellArray> (unstructuredgrid->GetCells());
      auto numCells = cellarray->GetNumberOfCells(); // vtkIdType
      this->mTriangles.reserve(numCells);
      // Traverse over VTK triangles
      cellarray->InitTraversal();
      auto idlist = vtkSmartPointer<vtkIdList>::New();
      while (cellarray->GetNextCell(idlist)) {
        if (idlist->GetNumberOfIds() != 3) // it is not a triangles
          continue;
        this->mTriangles.push_back({(size_t) idlist->GetId(0), (size_t) idlist->GetId(1), (size_t) idlist->GetId(2)});
      }
      // Shrink memory
      this->mPoints.shrink_to_fit();
      this->mTriangles.shrink_to_fit();
    }
  };
}}} // namespace

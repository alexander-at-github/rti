#pragma once

#include <vtkAbstractArray.h>
#include <vtkCellIterator.h>
#include <vtkCellData.h>
#include <vtkDataSetAttributes.h>
#include <vtkGenericCell.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtksys/SystemTools.hxx>
#include <vtkXMLPolyDataReader.h>

#include "rti/io/i_triangle_reader.hpp"

namespace rti { namespace io {
  template<typename Ty>
  class vtp_triangle_reader : public rti::io::i_triangle_reader<Ty> {
  public:
    vtp_triangle_reader(const std::string& pFilename) :
      rti::io::i_triangle_reader<Ty>(pFilename) {
      auto extension = vtksys::SystemTools::GetFilenameLastExtension(pFilename);
      auto extensionCondition = std::string(".vtp");
      if (extension != extensionCondition) {
        std::cerr
          << "Warning: " << typeid(this).name()
          << " may not be able to read the file " << pFilename
          << " because it does not have the " << extensionCondition
          << " extension." << std::endl;
      }
      auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
      reader->SetFileName(pFilename.c_str());
      reader->Update();
      vtkSmartPointer<vtkPolyData> polydata = reader->GetOutput();
      vtkIdType numPnts = polydata->GetNumberOfPoints();
      //vtkIdType numCells = polydata->GetNumberOfCells();

      auto pointdata = polydata->GetPointData();
      if (pointdata == nullptr) {
        std::cerr
          << "Warning: "  << typeid(this).name()
          << " could not find data in the file " << pFilename << std::endl;
      }
      this->mPoints.reserve(numPnts);
      for (vtkIdType idx = 0; idx < numPnts; ++idx) {
        double xyz[3]; // 3 dimensions
        polydata->GetPoint(idx, xyz);
        this->mPoints.push_back({(Ty) xyz[0], (Ty) xyz[1], (Ty) xyz[2]});
      }
      // Write triangles from the VTK data structure
      // for (vtkIdType idx = 0; idx < numCells; ++idx) {
      //   auto cell = polydata->GetCell(idx);
      // }
      auto cellarray = vtkSmartPointer<vtkCellArray> (polydata->GetPolys());
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
}} // namespace

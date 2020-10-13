#pragma once

#include <vtkAbstractArray.h>
#include <vtkCellIterator.h>
#include <vtkCellData.h>
#include <vtkDataSetAttributes.h>
#include <vtkGenericCell.h>
#include <vtkPointData.h>
//#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtksys/SystemTools.hxx>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>

#include "../i_point_cloud_reader.hpp"
#include "../../util/utils.hpp"

// There is an example at https://vtk.org/doc/nightly/html/classvtkCellIterator.html
//
// Good example at https://vtk.org/gitweb?p=VTK.git;a=blob;f=Examples/IO/Cxx/DumpXMLFile.cxx


namespace rti { namespace io { namespace xaver {
  // The parameter Ty is intended to be instantiated as a numeric type.
  template<typename Ty>
  class vtu_point_cloud_reader : public rti::io::i_point_cloud_reader<Ty> {
  public:
    vtu_point_cloud_reader(const std::string& pFilename) :
      mInfilename(pFilename) {
      std::string extension = vtksys::SystemTools::GetFilenameLastExtension(pFilename);
      if (extension != ".vtu") {
        std::cerr
          << "Warning: " << typeid(this).name()
          << " may not be able to read the file " << pFilename
          << " because it does not have the .vtp extension." << std::endl;
      }
      auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
      reader->SetFileName(pFilename.c_str());
      reader->Update();
      //vtkSmartPointer<vtkUnstructuredGrid> unstructuredgrid = reader->GetOutput(); // type: vtkUnstructuredGrid
      auto unstructuredgrid = reader->GetOutput(); // type: vtkUnstructuredGrid
      auto numPnts = unstructuredgrid->GetNumberOfPoints(); // type: vtkIdType
      auto pointdata = unstructuredgrid->GetPointData();
      auto celldata = unstructuredgrid->GetCellData();
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
      if (pointdata == nullptr || celldata == nullptr) {
        std::cerr
          << "Warning: "  << typeid(this).name()
          << " could not find data in the file " << pFilename << std::endl;
      }
      auto gsStr = "gridSpacing";
      auto gridSpacingArray = celldata->GetArray(gsStr);
      if (gridSpacingArray == nullptr) {
        std::cerr
          << "Warning: "  << typeid(this).name()
          << " could not find grid spacing data in the file " << pFilename << std::endl;
      }
      auto normalStr1 = "Normals";
      auto normalsArray = celldata->GetArray(normalStr1);
      if (normalsArray == nullptr) {
        auto normalStr2 = "normal_vector";
        normalsArray = pointdata->GetArray(normalStr2);
      }
      if (normalsArray == nullptr) {
        std::cerr
          << "Warning: "  << typeid(this).name()
          << " could not find surface normals in the file " << pFilename << std::endl;
      }
      for (vtkIdType idx = 0; idx < numPnts; ++idx) {
        double xyz[3]; // 3 dimensions
        unstructuredgrid->GetPoint(idx, xyz);
        double gridspacing[1]; // 1 dimensions
        gridSpacingArray->GetTuple(idx, gridspacing);
        auto radius = gridspacing[0];
        double radiusEpsilon = 1.0 / 4; // about 25%
        radius *= std::sqrt(3.0)/2 * (1 + radiusEpsilon);
        rti::util::quadruple<Ty> point {(Ty) xyz[0], (Ty) xyz[1] , (Ty) xyz[2], (Ty) radius};
        double nxnynz[3];
        normalsArray->GetTuple(idx, nxnynz);
        rti::util::triple<Ty> normal {(Ty) nxnynz[0], (Ty) nxnynz[1], (Ty) nxnynz[2]};
        // Do not add points which have radius equal to zero.
        // Do not add points which have a normal equal to the zero vector.
        if (radius == 0 || (nxnynz[0] == 0 && nxnynz[1] == 0 && nxnynz[2] == 0)) {
          continue;
        }
        mPoints.push_back(point);
        // Normalize
        if ( ! rti::util::is_normalized(normal))
          rti::util::normalize(normal);
        mNormals.push_back(normal);
      }
      // Shrink memory
      mPoints.shrink_to_fit();
      mNormals.shrink_to_fit();
    }

    std::vector<rti::util::quadruple<Ty> > get_points() override final {
      return mPoints;
    }

    std::vector<rti::util::triple<Ty> > get_normals() override final {
      return mNormals;
    }

    std::string get_input_file_name() const override final {
      return mInfilename;
    }
  private:
    std::string mInfilename;
    std::vector<rti::util::quadruple<Ty> > mPoints;
    std::vector<rti::util::triple<Ty> > mNormals;

  };
}}} // namespace

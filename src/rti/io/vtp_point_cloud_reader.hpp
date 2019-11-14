#pragma once

#include <boost/core/demangle.hpp>

#include <vtkAbstractArray.h>
#include <vtkCellIterator.h>
#include <vtkCellData.h>
#include <vtkDataSetAttributes.h>
#include <vtkGenericCell.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtksys/SystemTools.hxx>
#include <vtkXMLPolyDataReader.h>

#include "rti/io/i_point_cloud_reader.hpp"
#include "rti/util/utils.hpp"

// There is an example at https://vtk.org/doc/nightly/html/classvtkCellIterator.html
//
// Good example at https://vtk.org/gitweb?p=VTK.git;a=blob;f=Examples/IO/Cxx/DumpXMLFile.cxx


namespace rti { namespace io {
  // The parameter Ty is intended to be instantiated as a numeric type.
  template<typename Ty>
  class vtp_point_cloud_reader : public rti::io::i_point_cloud_reader<Ty> {
  public:
    vtp_point_cloud_reader(const std::string& pFilename) :
      mInfilename(pFilename) {
      std::string extension = vtksys::SystemTools::GetFilenameLastExtension(pFilename);
      if (extension != ".vtp") {
        std::cerr
          << "Warning: " << boost::core::demangle(typeid(this).name())
          << " may not be able to read the file " << pFilename
          << " because it does not have the .vtp extension." << std::endl;
      }
      auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
      reader->SetFileName(pFilename.c_str());
      reader->Update();
      vtkSmartPointer<vtkPolyData> polydata = reader->GetOutput();
      vtkIdType numPnts = polydata->GetNumberOfPoints();
      vtkIdType numCells = polydata->GetNumberOfCells();
      if (numPnts != numCells) {
        std::cerr
          << "Warning: " << boost::core::demangle(typeid(this).name())
          << "expects an input such that number of points equals number of cells. "
          << "File: " << pFilename << std::endl;
      }
      vtkSmartPointer<vtkCellData> celldata = polydata->GetCellData();
      // TODO: Clean up
      std::string sqrtOfAreaStr = "sqrt-of-area";
      celldata->SetActiveAttribute(sqrtOfAreaStr.c_str(), vtkDataSetAttributes::SCALARS);
      vtkSmartPointer<vtkDataArray> sqrtOfAreaArray = celldata->GetScalars();
      if (sqrtOfAreaArray == nullptr) {
        std::cerr
          << "Warning: " << boost::core::demangle(typeid(this).name())
          << " could not find cell data with the name \"sqrt-of-area\" in the file "
          << pFilename << std::endl;
      }
      vtkSmartPointer<vtkDataArray> normals = celldata->GetNormals();
      if (normals == nullptr) {
        std::cerr
          << "Warning: " << boost::core::demangle(typeid(this).name())
          << " could not find surface normal data in the file " << pFilename << std::endl;
      }
      for (vtkIdType idx = 0; idx < numPnts; ++idx) {
        double xyz[3]; // 3 dimensions
        polydata->GetPoint(idx, xyz);
        double radius[1]; // 1 dimension
        sqrtOfAreaArray->GetTuple(idx, radius);
        double radiusEpsilon = 1.0 / 32; // about 3%
        // TODO: Is that correct?
        radius[0] *= std::sqrt(3.0)/2 * (1 + radiusEpsilon);
        rti::util::quadruple<Ty> point {(Ty) xyz[0], (Ty) xyz[1] , (Ty) xyz[2], (Ty) radius[0]};
        mPoints.push_back(point);
        double nxnynz[3];
        normals->GetTuple(idx, nxnynz);
        rti::util::triple<Ty> normal {(Ty) nxnynz[0], (Ty) nxnynz[1], (Ty) nxnynz[2]};
        // Normalize
        if ( ! rti::util::is_normalized(normal))
          rti::util::normalize(normal);
        mNormals.push_back(normal);
      }
      // Shrink memory
      mPoints.shrink_to_fit();
      mNormals.shrink_to_fit();
    }

    ~vtp_point_cloud_reader() {}

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
}} // namespace rti

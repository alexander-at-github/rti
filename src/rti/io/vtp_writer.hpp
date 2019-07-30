#pragma once

#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedIntArray.h>
#include <vtkXMLPolyDataWriter.h>

#include "rti/geo/point_cloud_geometry.hpp"
#include "rti/geo/point_cloud_geometry.hpp"
#include "rti/trace/i_hit_counter.hpp"

namespace rti { namespace io {
  template<typename Ty>
  class vtp_writer {
  public:

    static
    bool write(rti::geo::point_cloud_geometry<Ty>& pGeometry,
               rti::trace::i_hit_counter& pHC,
               std::string pOutfilename) {
      // Precondition:
      assert (pGeometry.get_num_primitives() == pHC.get_counts().size() &&
              "hit count accumulator does not fit the given geometry");
      auto polydata = get_polydata(pGeometry);
      add_hit_counts(polydata, pHC);

      auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      writer->SetFileName(pOutfilename.c_str());
      writer->SetInputData(polydata);
      writer->SetDataModeToAscii(); // human readable XML output
      writer->Write();
      return true;
    }

  private:

    static
    void add_hit_counts(vtkSmartPointer<vtkPolyData> pPolydata,
                        rti::trace::i_hit_counter& pHC) {
      auto hitcnts = pHC.get_counts();
      assert (pPolydata->GetNumberOfPoints() == hitcnts.size() &&
              "polydata does not fit the hit count accumulator");
      auto hitValues = vtkSmartPointer<vtkUnsignedIntArray>::New();
      hitValues->SetNumberOfComponents(1); // 1 dimension
      hitValues->SetNumberOfTuples(hitcnts.size());
      for (size_t idx = 0; idx < hitcnts.size(); ++idx) {
        hitValues->InsertValue(idx, hitcnts[idx]);
      }
      hitValues->SetName(valueStr);
      pPolydata->GetCellData()->AddArray(hitValues);
    }

    static
    vtkSmartPointer<vtkPolyData> get_polydata(rti::geo::point_cloud_geometry<Ty>& pGeometry) {
      auto points = vtkSmartPointer<vtkPoints>::New();
      auto cells = vtkSmartPointer<vtkCellArray>::New();

      auto numpoints = pGeometry.get_num_primitives();
      //auto normals = std::vector<rti::util::triple<Ty> > {numpoints};
      auto normals = vtkSmartPointer<vtkDoubleArray>::New();
      normals->SetNumberOfComponents(3); // 3 dimensions
      normals->SetNumberOfTuples(numpoints);
      //auto radii = std::vector<Ty> {numpoints};
      auto radii = vtkSmartPointer<vtkDoubleArray>::New();
      radii->SetNumberOfComponents(1); // 1 dimensions
      radii->SetNumberOfTuples(numpoints);
      for (size_t idx = 0; idx < numpoints; ++idx) {
        auto point = pGeometry.get_prim(idx);
        auto normal = pGeometry.get_normal(idx);
        auto writePointId = points->InsertNextPoint(point[0], point[1], point[2]);
        cells->InsertNextCell(1, &writePointId); // one cell for writePointId
        //normals.push_back(normal);
        normals->SetTuple(idx, normal.data());
        //radii.push_back(point[3]);
        // the radius is saved as the 4th element of the point (which is an array)
        radii->SetTuple(idx, &point[3]);
      }
      auto polydata = vtkSmartPointer<vtkPolyData>::New();
      polydata->SetPoints(points);
      polydata->SetVerts(cells);
      polydata->GetCellData()->SetNormals(normals);
      radii->SetName(radiusStr);
      polydata->GetCellData()->AddArray(radii);

      return polydata;
    }
    
    // static
    // vtkSmartPointer<vtkPolyData> get_polydata(rti::geo::i_geometry<Ty>& pGeometry) {
    //   assert (false && "not implemented");
    //   return nullptr;
    // }

    static constexpr char const* radiusStr = "radius";
    static constexpr char const* valueStr = "value";
  };
}} // namespace

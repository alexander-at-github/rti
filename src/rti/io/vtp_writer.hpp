#pragma once

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkLine.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkXMLPolyDataWriter.h>

#include "rti/geo/i_boundary.hpp"
#include "rti/geo/point_cloud_sphere_geometry.hpp"
#include "rti/geo/point_cloud_disc_geometry.hpp"
#include "rti/trace/i_hit_accumulator.hpp"

namespace rti { namespace io {
  template<typename Ty>
  class vtp_writer {
  public:

    static
    void write(rti::geo::absc_point_cloud_geometry<Ty>& pGeometry,
               rti::trace::i_hit_accumulator<Ty>& pHC,
               std::string pOutfilename) {
      // Precondition:
      assert (pGeometry.get_num_primitives() == pHC.get_counts().size() &&
              "hit count accumulator does not fit the given geometry");
      auto polydata = get_polydata(pGeometry);
      add_hit_counts(polydata, pHC);
      write(polydata, pOutfilename);
    }

    static
    void write(rti::geo::i_boundary<Ty>& pBoundary,
               std::string pOutfilename) {
      auto polydata = get_polydata(pBoundary);
      write(polydata, pOutfilename);
    }

    static
    void write(std::vector<rti::util::pair<rti::util::triple<Ty> > >* pVec,
               std::string pOutfilename) {
      auto polydata = get_polydata(*pVec);
      write(polydata, pOutfilename);
    }

    static
    void write(std::vector<rti::util::triple<Ty> >* pVec,
               std::string pOutfilename) {
      auto polydata = get_polydata(*pVec);
      write(polydata, pOutfilename);
    }

  private:

    static
    void write(vtkSmartPointer<vtkPolyData> pPolydata, std::string pOutfilename) {
      auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      writer->SetFileName(pOutfilename.c_str());
      writer->SetInputData(pPolydata);
      writer->SetDataModeToAscii(); // human readable XML output
      writer->Write();
    }

    static
    void add_hit_counts(vtkSmartPointer<vtkPolyData> pPolydata,
                        rti::trace::i_hit_accumulator<Ty>& pAc) {
      auto hitcnts = pAc.get_counts();
      assert (pPolydata->GetNumberOfPoints() == hitcnts.size() &&
              "polydata does not fit the hit count accumulator");
      auto hitValues = vtkSmartPointer<vtkDoubleArray>::New();
      hitValues->SetNumberOfComponents(1); // 1 dimension
      hitValues->SetNumberOfTuples(hitcnts.size());
      for (size_t idx = 0; idx < hitcnts.size(); ++idx) {
        //if (hitcnts[idx] != 0) RLOG_DEBUG << "writing " << hitcnts[idx] << " to vtkPolyData" << std::endl;
        hitValues->InsertValue(idx, hitcnts[idx]);
      }
      hitValues->SetName(valueStr);
      pPolydata->GetCellData()->AddArray(hitValues);
    }

    static
    vtkSmartPointer<vtkPolyData> get_polydata(rti::geo::absc_point_cloud_geometry<Ty>& pGeometry) {
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

    static
    vtkSmartPointer<vtkPolyData> get_polydata(rti::geo::i_boundary<Ty>& pBoundary) {
      auto pointsOut = vtkSmartPointer<vtkPoints>::New();
      auto trianglesOut = vtkSmartPointer<vtkCellArray>::New();
      auto normalsOut = vtkSmartPointer<vtkDoubleArray>::New();
      for (auto const& pin : pBoundary.get_vertices()) {
        pointsOut->InsertNextPoint(pin[0], pin[1], pin[2]);
      }
      for (auto const& tin : pBoundary.get_triangles()) {
        auto triangleOut = vtkSmartPointer<vtkTriangle>::New();
        triangleOut->GetPointIds()->SetId(0, tin[0]);
        triangleOut->GetPointIds()->SetId(1, tin[1]);
        triangleOut->GetPointIds()->SetId(2, tin[2]);
        trianglesOut->InsertNextCell(triangleOut);
      }
      normalsOut->SetNumberOfComponents(3); // 3 dimensions
      normalsOut->SetNumberOfTuples(trianglesOut->GetNumberOfCells());
      // for (auto const& nin : pBoundary.get_triangle_normals()) {
      //   normalsOut->InsertNextTuple(nin.data()); // 0 is the start index of the second argument
      // }
      auto nin = pBoundary.get_triangle_normals();
      for (size_t idx = 0; idx < nin.size(); ++idx) {
        normalsOut->SetTuple(idx, nin[idx].data());
      }

      auto polydata = vtkSmartPointer<vtkPolyData>::New();
      polydata->SetPoints(pointsOut);
      polydata->SetPolys(trianglesOut);
      polydata->GetCellData()->SetNormals(normalsOut);
      return polydata;
    }

    static
    vtkSmartPointer<vtkPolyData> get_polydata(
        std::vector<rti::util::pair<rti::util::triple<float> > >& pVec) {
      auto pointsOut = vtkSmartPointer<vtkPoints>::New();
      auto linesOut = vtkSmartPointer<vtkCellArray>::New();

      auto pidx = 0ull;
      for (auto const& lineIn : pVec) {
        pointsOut->InsertNextPoint(lineIn[0].data());
        pointsOut->InsertNextPoint(lineIn[1].data());
        auto lineOut = vtkSmartPointer<vtkLine>::New();
        lineOut->GetPointIds()->SetId(0,pidx+0);
        lineOut->GetPointIds()->SetId(1,pidx+1);
        linesOut->InsertNextCell(lineOut);
        pidx += 2; // 2 points have been inserted
      }
      assert (pidx == 2 * pVec.size() && "Missmatch between point-index and size of vector");

      auto polydata = vtkSmartPointer<vtkPolyData>::New();
      polydata->SetPoints(pointsOut);
      polydata->SetLines(linesOut);
      return polydata;
    }

    static
    vtkSmartPointer<vtkPolyData> get_polydata(
        std::vector<rti::util::triple<float> >& pVec) {
      auto pointsOut =vtkSmartPointer<vtkPoints>::New();
      for (auto const& pointIn : pVec) {
        pointsOut->InsertNextPoint(pointIn.data());
      }
      auto polydata = vtkSmartPointer<vtkPolyData>::New();
      polydata->SetPoints(pointsOut);
      return polydata;
    }

    static constexpr char const* radiusStr = "radius";
    static constexpr char const* valueStr = "value";
  };
}} // namespace

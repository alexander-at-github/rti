#include <boost/core/demangle.hpp> // debug // remove!
#include <boost/algorithm/string.hpp>

#include <cmath>
#include <sstream>
#include <unordered_map>

#include "vtkActor.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkGlyph3D.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkTypeInt32Array.h"
//#include "vtkVertexGlyphFilter.h"
#include "vtkXMLPolyDataWriter.h"

#include "rti/enum_class_hash_function.hpp"

// An example how to write (or read) normals in vtk is at the following URL.
// https://vtk.org/Wiki/VTK/Examples/Cxx/PolyData/PolyDataCellNormals
// ( https://vtk.org/Wiki/VTK/Tutorials/TriangleGeometryVertices )

namespace rti::dsv {
  enum class cml_opt_type {INFILE_NAME, OUTFILE_NAME, FILTER_COVERED, DISPLAY};
  class command_line_options {
    // A simple command line options implementation
  public:
    command_line_options(int argc, char** argv) {
      // Start from one. argv[0] is equal to the file name of the executable.
      for (size_t idx = 1; idx < argc-1; ++idx) {
        for (auto& oo : options) {
          if (oo.second.optStr == argv[idx]) {
            oo.second.value = argv[idx+1];
          }
        }
      }
    }
    std::string get_value(rti::dsv::cml_opt_type pType) {
      for (auto& oo : options) {
        if (oo.first /* key of map */ == pType) {
          return oo.second.value;
        }
      }
      return "";
    }
  private:
    struct option_spec {
      std::string name;
      std::string optStr;
      std::string value;
    };
    std::unordered_map<cml_opt_type, option_spec, rti::enum_class_hash_function> options {
      {cml_opt_type::INFILE_NAME, option_spec {"input file name", "--infile", ""}},
      {cml_opt_type::OUTFILE_NAME, option_spec {"output file name", "--write", ""}},
      {cml_opt_type::FILTER_COVERED, option_spec {"filter covered points", "--filter-covered", "false"}},
      {cml_opt_type::DISPLAY, option_spec{"display the input on GUI", "--display", "false"}}
    };
  };

} // namespace rti::dsv

void print_usage(std::string pName) {
  std::cout
    << "Usage: " << pName
    << " [--write <outfile-name>] [--filter-covered [true | false]] [--display [true | false]] --infile <infile-name>" << std::endl;
}


// Partially built on the ReadTextFile example in VTK
int main(int argc, char* argv[]) {
  rti::dsv::command_line_options options(argc, argv);
  if (options.get_value(rti::dsv::cml_opt_type::INFILE_NAME) == "") {
    // No input file given
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }

  std::string infilename = options.get_value(rti::dsv::cml_opt_type::INFILE_NAME);
  std::string outfilename = options.get_value(rti::dsv::cml_opt_type::OUTFILE_NAME);
  bool filterCovered =
    options.get_value(rti::dsv::cml_opt_type::FILTER_COVERED) == "true" ? true : false;
  bool display =
    options.get_value(rti::dsv::cml_opt_type::DISPLAY) == "true" ? true : false;
  if (filterCovered) {
    std::cout << "Filtering points with cover flag set to a value not equal zero." << std::endl;
  }

  std::ifstream filestream(infilename.c_str());

  // Read
  std::string line;
  // VTK calls topological elements points and geometric (visual) elements vertices
  auto points = vtkSmartPointer<vtkPoints>::New();
  // Array of vertex cells
  auto vertices = vtkSmartPointer<vtkCellArray>::New();
  // Temporary storage for surface normals and other fields in input file
  std::vector<std::array<double, 3> > normalList;
  std::vector<int32_t> matIdList;
  std::vector<double> areaList;
  std::vector<int32_t> coverList;
  while (std::getline(filestream, line)) {
    boost::trim(line);
    if (line[0] == '#') {
      // The string in the variable line represents a comment
      continue;
    }
    std::stringstream linestream;
    linestream << line;
    double xx, yy, zz;
    linestream >> xx >> yy >> zz;
    double nx, ny, nz;
    int32_t mid;
    double area;
    int32_t cover;
    linestream >> nx >> ny >> nz >> mid >> area >> cover;
    if (filterCovered && cover != 0) {
      // Skip that point. It is covered by another point on a finer level
      continue;
    }
    auto pointId = points->InsertNextPoint(xx, yy, zz);
    // create a vertex cell on the point
    vertices->InsertNextCell(1, &pointId); // arguments: number of cells, array of point ids
    normalList.push_back({nx, ny, nz});
    matIdList.push_back(mid);
    areaList.push_back(area);
    coverList.push_back(cover);
  }
  filestream.close();
  auto polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);
  polydata->SetVerts(vertices);


  // Set vertex (cell) normals
  auto normalsVtkArray = vtkSmartPointer<vtkDoubleArray>::New();
  normalsVtkArray->SetNumberOfComponents(3); // 3 dimensions
  normalsVtkArray->SetNumberOfTuples(normalList.size());
  for (size_t idx = 0; idx < normalList.size(); ++idx) {
    normalsVtkArray->SetTuple(idx, normalList[idx].data());
  }
  // add normals to the cells in the polydata
  polydata->GetCellData()->SetNormals(normalsVtkArray);

  // Set other data
  auto matIdVtkArray = vtkSmartPointer<vtkTypeInt32Array>::New();
  auto areaVtkArray = vtkSmartPointer<vtkDoubleArray>::New();
  auto sqrtOfAreaVtkArray = vtkSmartPointer<vtkDoubleArray>::New();
  auto coverVtkArray = vtkSmartPointer<vtkTypeInt32Array>::New();
  matIdVtkArray->SetName("mat-id");
  areaVtkArray->SetName("area");
  sqrtOfAreaVtkArray->SetName("sqrt-of-area");
  coverVtkArray->SetName("cover-flag");
  matIdVtkArray->SetNumberOfComponents(1);
  areaVtkArray->SetNumberOfComponents(1);
  sqrtOfAreaVtkArray->SetNumberOfComponents(1);
  coverVtkArray->SetNumberOfComponents(1);
  matIdVtkArray->SetNumberOfTuples(matIdList.size());
  areaVtkArray->SetNumberOfTuples(matIdList.size());
  sqrtOfAreaVtkArray->SetNumberOfTuples(matIdList.size());
  coverVtkArray->SetNumberOfTuples(matIdList.size());
  for (size_t idx = 0; idx < matIdList.size(); ++idx) {
    matIdVtkArray->SetTypedTuple(idx, &matIdList[idx]);
    areaVtkArray->SetTuple(idx, &areaList[idx]);
    double sqrtOfArea = std::sqrt(areaList[idx]);
    sqrtOfAreaVtkArray->SetTuple(idx, &sqrtOfArea);
    coverVtkArray->SetTypedTuple(idx, &coverList[idx]);
  }
  polydata->GetCellData()->AddArray(matIdVtkArray);
  polydata->GetCellData()->AddArray(areaVtkArray);
  polydata->GetCellData()->AddArray(sqrtOfAreaVtkArray);
  polydata->GetCellData()->AddArray(coverVtkArray);

  if (display) {
    // Visualize
    //auto glyphFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    //glyphFilter->SetInputData(polydata);
    //glyphFilter->Update();

    auto sphere = vtkSmartPointer<vtkSphereSource>::New();
    //sphere->SetRadius(0.07); // magic number
    sphere->SetRadius(0.8); // magic number
    auto glyph = vtkSmartPointer<vtkGlyph3D>::New();
    //glyph->SetInputConnection(glyphFilter->GetOutputPort());
    glyph->SetInputData(polydata);
    glyph->Update();
    glyph->SetSourceConnection(sphere->GetOutputPort());
    glyph->ScalingOn(); // necessary? NO
    glyph->SetScaleModeToScaleByScalar();
    //glyph->SetColorModeToColorByScale();
    glyph->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "sqrt-of-area");

    //
    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(glyph->GetOutputPort());
    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    auto renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->AddActor(actor);
    renderer->SetBackground(0.3, 0.6, 0.3);
    auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    auto renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);
    renderWindowInteractor->Start();
  }

  if (outfilename != "") {
    // Write file
    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(outfilename.c_str());
    writer->SetInputData(polydata);
    writer->SetDataModeToAscii(); // human readable XML output
    writer->Write();
    std::cout << "File " << outfilename << " written." << std::endl;
  }

  return EXIT_SUCCESS;
}

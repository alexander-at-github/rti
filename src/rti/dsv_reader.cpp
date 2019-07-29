#include <boost/algorithm/string.hpp>

#include <cmath>
#include <sstream>
#include <unordered_map>

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkGlyph3D.h>
#include <vtkPoints.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTypeInt32Array.h>
//#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPolyDataWriter.h>

#include "rti/util/clo.hpp"
#include "rti/util/enum_class_hash_function.hpp"

// An example how to write (or read) normals in vtk is at the following URL.
// https://vtk.org/Wiki/VTK/Examples/Cxx/PolyData/PolyDataCellNormals
// ( https://vtk.org/Wiki/VTK/Tutorials/TriangleGeometryVertices )

// Partially built on the ReadTextFile example in VTK
int main(int argc, char* argv[]) {

  auto optMan = std::make_unique<rti::util::clo::manager>();
  optMan->addCmlParam(rti::util::clo::bool_option
    {"FILTER_COVERED", {"--filter-covered"}, "turns filtering of covered points on"});
  optMan->addCmlParam(rti::util::clo::bool_option
    {"RENDER", {"--render"}, "render the input on GUI"});
  optMan->addCmlParam(rti::util::clo::string_option
    {"INPUT_FILE", {"--infile"}, "spacifies the name of the input file", true});
  optMan->addCmlParam(rti::util::clo::string_option
    {"OUTPUT_FILE", {"--write", "--output"}, "specifies the name of the output file", false});
  bool succ = optMan->parse_args(argc, argv);
  if (!succ) {
    std::cout << optMan->get_usage_msg();
    exit(EXIT_FAILURE);
  }
  std::string infilename = optMan->get_string_option_value("INPUT_FILE");
  std::string outfilename = optMan->get_string_option_value("OUTPUT_FILE");
  bool filterCovered = optMan->get_bool_option_value("FILTER_COVERED");
  bool render = optMan->get_bool_option_value("RENDER");

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
    double nx, ny, nz;
    int32_t mid;
    double area;
    int32_t cover;
    linestream >> xx >> yy >> zz >> nx >> ny >> nz >> mid >> area >> cover;
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
  // One could also do polydata->GetPointData()->AddArray(...)
  polydata->GetCellData()->AddArray(matIdVtkArray);
  polydata->GetCellData()->AddArray(areaVtkArray);
  polydata->GetCellData()->AddArray(sqrtOfAreaVtkArray);
  polydata->GetCellData()->AddArray(coverVtkArray);

  if (render) {
    // Visualize
    //auto glyphFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    //glyphFilter->SetInputData(polydata);
    //glyphFilter->Update();

    auto sphere = vtkSmartPointer<vtkSphereSource>::New();
    sphere->SetThetaResolution(16); // magic number; only for rendering
    sphere->SetPhiResolution(16); // magic number; only for rendering
    //sphere->SetRadius(0.07); // magic number
    sphere->SetRadius(std::sqrt(3) / 2 * (1 + 1/32 /* increase by about 3% */)); // magic number
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

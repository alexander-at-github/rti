// Based on the ReadTextFile example in VTK

#include <boost/algorithm/string.hpp>

#include <sstream>

#include "vtkActor.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSmartPointer.h"
#include "vtkTypeInt32Array.h"
#include "vtkVertexGlyphFilter.h"
#include "vtkXMLPolyDataWriter.h"

// An example how to write (or read) normals in vtk is at the following URL.
// https://vtk.org/Wiki/VTK/Examples/Cxx/PolyData/PolyDataCellNormals
// ( https://vtk.org/Wiki/VTK/Tutorials/TriangleGeometryVertices )


void print_usage(std::string pName) {
  std::cout << "Usage: " << pName << " [--write <outfile-name>] <infile-name>" << std::endl;
}

int main(int argc, char* argv[]) {
  std::cerr << "argc: " << argc << std::endl;
  if ( ! (argc == 2 || (argc == 4 && std::string(argv[1]) == "--write"))) {
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }
  std::string outfilename;
  if (argc == 4) {
    outfilename = argv[2];
  }

  std::string infilename = argv[argc-1];
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
    double xx, yy, zz;
    std::stringstream linestream;
    linestream << line;
    linestream >> xx >> yy >> zz;
    auto pointId = points->InsertNextPoint(xx, yy, zz);
    // create a vertex cell on the point
    vertices->InsertNextCell(1, &pointId); // arguments: number of cells, array of point ids
    double nx, ny, nz;
    int32_t mid;
    double area;
    int32_t cover;
    linestream >> nx >> ny >> nz >> mid >> area >> cover;
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
  auto coverVtkArray = vtkSmartPointer<vtkTypeInt32Array>::New();
  matIdVtkArray->SetName("mat-id");
  areaVtkArray->SetName("area");
  coverVtkArray->SetName("cover-flag");
  matIdVtkArray->SetNumberOfComponents(1);
  areaVtkArray->SetNumberOfComponents(1);
  coverVtkArray->SetNumberOfComponents(1);
  matIdVtkArray->SetNumberOfTuples(matIdList.size());
  areaVtkArray->SetNumberOfTuples(matIdList.size());
  coverVtkArray->SetNumberOfTuples(matIdList.size());
  for (size_t idx = 0; idx < matIdList.size(); ++idx) {
    matIdVtkArray->SetTypedTuple(idx, &matIdList[idx]);
    areaVtkArray->SetTuple(idx, &areaList[idx]);
    coverVtkArray->SetTypedTuple(idx, &coverList[idx]);
  }
  polydata->GetCellData()->AddArray(matIdVtkArray);
  polydata->GetCellData()->AddArray(areaVtkArray);
  polydata->GetCellData()->AddArray(coverVtkArray);

  auto glyphFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
  glyphFilter->SetInputData(polydata);
  glyphFilter->Update();

  // Visualize
  auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(glyphFilter->GetOutputPort());
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

// Based on the ReadTextFile example in VTK

#include <sstream>

#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkVertexGlyphFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"


int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " Filename" << std::endl;
    return EXIT_FAILURE;
  }
  std::string filename = argv[1];
  std::ifstream filestream(filename.c_str());

  // Read
  std::string line;
  auto points = vtkSmartPointer<vtkPoints>::New();
  while (std::getline(filestream, line)) {
    double xx, yy, zz;
    std::stringstream linestream;
    linestream << line;
    linestream >> xx >> yy >> zz;
    points->InsertNextPoint(xx, yy, zz);
  }
  filestream.close();
  auto polyData = vtkSmartPointer<vtkPolyData>::New();
  polyData->SetPoints(points);
  auto glyphFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
  glyphFilter->SetInputData(polyData);
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


  return EXIT_SUCCESS;
}

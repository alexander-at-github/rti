#include "rti/io/gmsh_reader.hpp"
#include "rti/util/clo.hpp"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>

namespace rti { namespace main_gmsh_reader {
  rti::util::triple<double> get_normal(int const& pPointid,
                                       std::vector<rti::util::triple<double> >& pPoints,
                                       std::vector<rti::util::triple<std::size_t> >& pTriangles) {
    auto normals = std::vector<rti::util::triple<double> > {};
    for (auto& triangle : pTriangles) {
      if (rti::util::contains<std::size_t>(triangle, pPointid)) {
        auto trinum = rti::util::triple<rti::util::triple<double> > {
          pPoints[triangle[0]], pPoints[triangle[1]], pPoints[triangle[2]]};
        normals.push_back(rti::util::compute_normal<double>(trinum));
      }
    }
    auto result = rti::util::triple<double> {0, 0, 0};
    for (auto const& nn : normals) {
      result[0] += nn[0];
      result[1] += nn[1];
      result[2] += nn[2];
    }
    rti::util::normalize<double>(result);
    return result;
  }

    double get_radius(int const& pPointid,
                      std::vector<rti::util::triple<double> >& pPoints,
                      std::vector<rti::util::triple<std::size_t> >& pTriangles) {
      auto radius = 0.0; // double
      for (auto& triangle : pTriangles) {
        if (rti::util::contains<std::size_t>(triangle, pPointid)) {
          auto centroid = rti::util::centroid<double>({pPoints[triangle[0]],
                                                       pPoints[triangle[1]],
                                                       pPoints[triangle[2]]});
          auto distance = rti::util::distance<double>({pPoints[pPointid], centroid});
          radius = std::max<double> (radius, distance);
        }
      }
      return radius;
    }
}} // namespace


int main(int argc, char* argv[]) {
  auto optman = rti::util::clo::manager {};
  optman.addCmlParam(rti::util::clo::string_option
                     {"INPUT_FILE", {"--infile"},
                        "spacifies the name of the input file", true});
  optman.addCmlParam(rti::util::clo::string_option
                     {"OUTPUT_FILE", {"--outfile", "--write"},
                        "specifies the name of the output file", true});
  optman.addCmlParam(rti::util::clo::bool_option
                     {"INVERT_NORMALS", {"--invert-normals"},
                        "inverts the surface normals"});
  optman.addCmlParam(rti::util::clo::string_option
                     {"LOWER_BOUND", {"--lower-bound"},
                        "remove all points with an z-value less than this argument", false});
  optman.addCmlParam(rti::util::clo::string_option
                     {"UPPER_BOUND", {"--upper-bound"},
                        "remove all points with an z-value higher than this argument", false});
  auto succ = optman.parse_args(argc, argv);
  if (!succ) {
    std::cout << optman.get_usage_msg();
    exit(EXIT_FAILURE);
  }
  auto infilename = optman.get_string_option_value("INPUT_FILE");
  auto outfilename = optman.get_string_option_value("OUTPUT_FILE");
  auto invertnormals = optman.get_bool_option_value("INVERT_NORMALS");
  auto lowerboundstr = optman.get_string_option_value("LOWER_BOUND");
  auto lowerbound = lowerboundstr.empty() ? std::numeric_limits<double>::lowest() : std::stod(lowerboundstr);
  auto upperboundstr = optman.get_string_option_value("UPPER_BOUND");
  auto upperbound = upperboundstr.empty() ? std::numeric_limits<double>::max() : std::stod(upperboundstr);
  // std::cerr << "lowerboundstr == " << lowerboundstr << std::endl;
  // std::cerr << "upperboundstr == " << upperboundstr << std::endl;
  // std::cerr << "lowerbound == " << lowerbound << std::endl;
  // std::cerr << "upperbound == " << upperbound << std::endl;

  auto& reader = rti::io::gmsh_reader::getInstance(infilename);

  // VTK calls topological elements points and geometric (visual) elements vertices
  auto points = vtkSmartPointer<vtkPoints>::New();
  // Array of vertex cells
  auto vertices = vtkSmartPointer<vtkCellArray>::New();

  // Read
  auto inputpoints = reader.get_vertices();
  auto inputtriangles = reader.get_triangles();
  std::cerr << "inputpoints.size() == " << inputpoints.size() << std::endl;
  std::cerr << "lowerbound == " << lowerbound << std::endl;
  std::cerr << "upperbound == " << upperbound << std::endl;

  // Process input
  auto normallist = std::vector<rti::util::triple<double> > {};
  auto radiuslist = std::vector<double> {};
  for (size_t idx = 0; idx < inputpoints.size(); ++idx) {
    //for (auto const& inpoint : inputpoints) {
    auto const& inpoint = inputpoints[idx];
    //std::cerr << "inpoint == " << inpoint[0] << " " << inpoint[1] << " " << inpoint[2] << std::endl;
    if ( ! (lowerbound <= inpoint[2] && inpoint[2] <= upperbound)) { // filter z-axis
      continue;
    }
    auto normal = rti::main_gmsh_reader::get_normal(idx, inputpoints, inputtriangles);
    auto radius = rti::main_gmsh_reader::get_radius(idx, inputpoints, inputtriangles) * 1.5; /* plus 10 percent */
    if (invertnormals)
      normal = rti::util::inv<double>(normal);
    auto pointid = points->InsertNextPoint(inpoint[0], inpoint[1], inpoint[2]);
    vertices->InsertNextCell(1, &pointid);
    normallist.push_back({normal[0], normal[1], normal[2]});
    radiuslist.push_back(radius);
  }
  std::cerr << "normallist.size() == " << normallist.size() << std::endl;
  auto polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);
  polydata->SetVerts(vertices);
  // Set Normals
  auto normalsVtkArray = vtkSmartPointer<vtkDoubleArray>::New();
  normalsVtkArray->SetNumberOfComponents(3); // 3 dimensions
  normalsVtkArray->SetNumberOfTuples(normallist.size());
  for (size_t idx = 0; idx < normallist.size(); ++idx) {
    normalsVtkArray->SetTuple(idx, normallist[idx].data());
  }
  polydata->GetCellData()->SetNormals(normalsVtkArray);
  // Set area / radius / square root of area
  auto sqrtOfAreaVtkArray = vtkSmartPointer<vtkDoubleArray>::New();
  sqrtOfAreaVtkArray->SetName("sqrt-of-area");
  sqrtOfAreaVtkArray->SetNumberOfComponents(1);
  sqrtOfAreaVtkArray->SetNumberOfTuples(radiuslist.size());
  for (size_t idx = 0; idx < radiuslist.size(); ++idx)
    sqrtOfAreaVtkArray->SetTuple(idx, &radiuslist[idx]);
  polydata->GetCellData()->AddArray(sqrtOfAreaVtkArray);

  // Write to file
  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(outfilename.c_str());
  writer->SetInputData(polydata);
  writer->SetDataModeToAscii(); // human readable XML output
  writer->Write();
  if (invertnormals)
    std::cout << "Inverted surface normals" << std::endl;
  std::cout << "File " << outfilename << " written." << std::endl;

  return EXIT_SUCCESS;
}

#include <config.h>

#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/curvedgeometry/geometries/sphere.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/gmsh4/gmsh4reader.hh>

#include "mean_curvature_flow.hh"
#include "runner.hh"

int main(int argc, char *argv[])
{
  using namespace Dune;
  MPIHelper::instance(argc, argv);

  std::string inifile = "mcf.ini";
  if (argc > 1)
    inifile = argv[1];

  ParameterTree pt;
  ParameterTreeParser::readINITree(inifile, pt);

  std::string gridFilename = pt.get<std::string>("grid.filename", DUNE_GRID_PATH "sphere_very_rough.msh")

  using HostGrid = FoamGrid<2,3>;
  auto hostGridPtr = Gmsh4Reader<HostGrid>::createGridFromFile(gridFilename);

  double radius = pt.get<double>("grid.initial.radius", 1.0);
  auto initialSurface = sphereGridFunction(radius);

  int kg = pt.get<int>("grid.kg", 2);
  switch (kg) {
  case 1: run<1>(pt, *hostGridPtr, initialSurface, BGN::MeanCurvatureFlow{}); break;
  case 2: run<2>(pt, *hostGridPtr, initialSurface, BGN::MeanCurvatureFlow{}); break;
  case 3: run<3>(pt, *hostGridPtr, initialSurface, BGN::MeanCurvatureFlow{}); break;
  default:
    DUNE_THROW(NotImplemented, "call run<kg>(...) for your polynomial order.");
  }
}

#include <config.h>

#include <cmath>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/curvedgeometry/geometries/sphere.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/geometry/types.hh>
#include <dune/gmsh4/gmsh4reader.hh>
#include <dune/grid/common/gridfactory.hh>

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

  using HostGrid = FoamGrid<1,2>;
  using Factory = GridFactory<HostGrid>;

  int refinement = pt.get<int>("grid.initial.refinement", 10);
  double radius = pt.get<double>("grid.initial.radius", 1.0);

  Factory factory;
  for (int i = 0; i < refinement; ++i) {
    double theta = i*2.0*M_PI/refinement;
    factory.insertVertex(FieldVector<double,2>{radius * std::cos(theta), radius * std::sin(theta)});
  }
  for (int i = 0; i < refinement; ++i) {
    factory.insertElement({i, (i+1)%refinement}, GeometryTypes::line);
  }
  auto hostGridPtr = factory.createGrid();

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

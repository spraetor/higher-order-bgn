#include <config.h>

#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/curvedgrid/geometries/sphere.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/gmsh4/gmsh4reader.hh>

#include <dune/higher-order-bgn/runner.hh>
#include <dune/higher-order-bgn/surfaceprojection.hh>

#if FLOW == 1
#include <dune/higher-order-bgn/meancurvatureflow.hh>
using GeometricFlow = Dune::BGN::MeanCurvatureFlow;
#elif FLOW == 2
#include <dune/higher-order-bgn/surfacediffusion.hh>
using GeometricFlow = Dune::BGN::SurfaceDiffusion;
#endif

int main(int argc, char *argv[])
{
  using namespace Dune;
  MPIHelper::instance(argc, argv);

  std::string inifile = "mcf.ini";
  if (argc > 1)
    inifile = argv[1];

  ParameterTree pt;
  ParameterTreeParser::readINITree(inifile, pt);

  std::string gridFilename = pt.get<std::string>("grid.filename", DUNE_GRID_PATH "sphere_very_rough.msh");

  using HostGrid = FoamGrid<2,3>;
  auto hostGridPtr = Gmsh4Reader<HostGrid>::createGridFromFile(gridFilename);

#if SURFACE == 1
  double radius = pt.get<double>("grid.initial.radius", 1.0);
  auto initialSurface = SphereProjection<3,double>{radius};
#elif SURFACE == 2
  double a = pt.get<double>("grid.initial.a", 1.0);
  double b = pt.get<double>("grid.initial.b", 1.0);
  double c = pt.get<double>("grid.initial.c", 1.0);
  auto initialSurface = BGN::EllipsoidProjection<double>{a,b,c};
#endif

  int kg = pt.get<int>("grid.kg", 2);
  switch (kg) {
  case 1: run<1>(pt, *hostGridPtr, initialSurface, GeometricFlow{}); break;
  case 2: run<2>(pt, *hostGridPtr, initialSurface, GeometricFlow{}); break;
  case 3: run<3>(pt, *hostGridPtr, initialSurface, GeometricFlow{}); break;
  default:
    DUNE_THROW(NotImplemented, "call eoc<kg>(...) for your polynomial order.");
  }
}

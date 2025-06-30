#include <config.h>

#include <cmath>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/curvedgrid/geometries/sphere.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/geometry/type.hh>
#include <dune/gmsh4/gmsh4reader.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dune/higher-order-bgn/curveprojection.hh>
#include <dune/higher-order-bgn/eoc.hh>

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

  using HostGrid = FoamGrid<1,2>;

  unsigned int refinement = pt.get<int>("grid.initial.refinement", 10);
#if SURFACE == 1
  double radius = pt.get<double>("grid.initial.radius", 1.0);
#elif SURFACE == 2
  double a = pt.get<double>("grid.initial.a", 1.0);
  double b = pt.get<double>("grid.initial.b", 1.0);
#endif

  // create a circle grid by explicitly inserting vertices and connectivity
  GridFactory<HostGrid> factory;
  for (unsigned int i = 0; i < refinement; ++i) {
    double theta = i*2.0*M_PI/refinement;
#if SURFACE == 1
    factory.insertVertex(FieldVector<double,2>{radius * std::cos(theta), radius * std::sin(theta)});
#elif SURFACE == 2
    factory.insertVertex(FieldVector<double,2>{a * std::cos(theta), b * std::sin(theta)});
#endif
  }
  for (unsigned int i = 0; i < refinement; ++i) {
    factory.insertElement(GeometryTypes::line, {i, (i+1)%refinement});
  }
  auto hostGridPtr = factory.createGrid();

#if SURFACE == 1
  auto initialSurface = SphereProjection<2,double>{radius};
#elif SURFACE == 2
  auto initialSurface = BGN::EllipseProjection<double>{a,b};
#endif

  int kg = pt.get<int>("grid.kg", 2);
  switch (kg) {
  case 1: eoc<1>(pt, *hostGridPtr, initialSurface, GeometricFlow{}); break;
  case 2: eoc<2>(pt, *hostGridPtr, initialSurface, GeometricFlow{}); break;
  case 3: eoc<3>(pt, *hostGridPtr, initialSurface, GeometricFlow{}); break;
  default:
    DUNE_THROW(NotImplemented, "call eoc<kg>(...) for your polynomial order.");
  }
}

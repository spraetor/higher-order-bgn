#include <config.h>

#include <cassert>
#include <cmath>
#include <thread>

#include <dune/common/indices.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/higher-order-bgn/boundarygridfactory.hh>
#include <dune/higher-order-bgn/eoc.hh>

#if FLOW == 1
#include <dune/higher-order-bgn/meancurvatureflow.hh>
using GeometricFlow = Dune::BGN::MeanCurvatureFlow;
#elif FLOW == 2
#include <dune/higher-order-bgn/surfacediffusion.hh>
using GeometricFlow = Dune::BGN::SurfaceDiffusion;
#endif

// Create a surface grid representing the boundary of an elongated box
template <class GridType>
auto createBoxBoundaryGrid(double dx, double dy, double dz)
{
  double d = std::min({dx,dy,dz});
  using BoxGridType = Dune::UGGrid<GridType::dimensionworld>;
  auto boxGrid = Dune::StructuredGridFactory<BoxGridType>::createSimplexGrid({0.0,0.0,0.0}, {dx,dy,dz}, {(unsigned int)(dx/d),(unsigned int)(dy/d),(unsigned int)(dz/d)});
  return Dune::BoundaryGridFactory<GridType>::createBoundaryGrid(boxGrid->leafGridView());
}

int main(int argc, char *argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  std::string inifile = "mcf.ini";
  if (argc > 1)
    inifile = argv[1];

  Dune::ParameterTree pt;
  Dune::ParameterTreeParser::readINITree(inifile, pt);

  auto hostGridPtr = createBoxBoundaryGrid<Dune::FoamGrid<2,3>>(10,1,1);
  auto identityMap = [](auto const& x) { return x; };

  int kg = pt.get<int>("grid.kg", 2);
  switch (kg) {
  case 1: BGN::eoc<1>(pt, *hostGridPtr, identityMap, GeometricFlow{}); break;
  case 2: BGN::eoc<2>(pt, *hostGridPtr, identityMap, GeometricFlow{}); break;
  case 3: BGN::eoc<3>(pt, *hostGridPtr, identityMap, GeometricFlow{}); break;
  default:
    DUNE_THROW(NotImplemented, "call eoc<kg>(...) for your polynomial order.");
  }
}

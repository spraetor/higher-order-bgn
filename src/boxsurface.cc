#include <config.h>

#include <cassert>
#include <cmath>
#include <thread>

#include <dune/common/indices.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/meshdist/hausdorffdistance.hh>

#include "boundarygridfactory.hh"
#include "gridsize.hh"
#include "mean_curvature_flow.hh"
#include "printerrors.hh"
#include "runner.hh"
#include "surface_diffusion.hh"

using namespace Dune;
using namespace std;

// Define the flow type as a local-assembler factory
using Flow = Dune::BGN::SurfaceDiffusion;
// using Float = MeanCurvatureFlow;

// Create a surface grid representing the boundary of an elongated box
template <class GridType>
auto createBoxBoundaryGrid(double dx, double dy, double dz)
{
  double d = std::min({dx,dy,dz});
  using BoxGridType = Dune::UGGrid<GridType::dimensionworld>;
  auto boxGrid = Dune::StructuredGridFactory<BoxGridType>::createSimplexGrid({0.0,0.0,0.0}, {dx,dy,dz}, {(unsigned int)(dx/d),(unsigned int)(dy/d),(unsigned int)(dz/d)});
  return Dune::BoundaryGridFactory<GridType>::createBoundaryGrid(boxGrid->leafGridView());
}

// Create a function-space basis [V^3 x V] suitable for the geometric flow problem
template <int k = 2, class GridView>
auto makeFlowBasis(GridView const& gridView)
{
  using namespace Dune::Functions::BasisFactory;
  return makeBasis(gridView,
      composite(power<GridView::dimensionworld>(lagrange<k>(), flatInterleaved()),
                lagrange<k>(),
                flatLexicographic()));
}

int main(int argc, char *argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  std::string inifile = "mcf.ini";
  if (argc > 1)
    inifile = argv[1];

  std::size_t threadCount = std::thread::hardware_concurrency();
  if (argc>2)
    threadCount = std::stoul(std::string(argv[2]));
  std::cout << "Using parallel executor with " << threadCount << " threads" << std::endl;

  Dune::ParameterTree pt;
  Dune::ParameterTreeParser::readINITree(inifile, pt);

  auto hostGrid = createBoxBoundaryGrid<Dune::FoamGrid<2,3>>(10,1,1);

  int maxLevel = pt.get<int>("grid.refinement_levels", 3);
  hostGrid->globalRefine(maxLevel+2);

  auto identityMap = [](auto const& x) { return x; };

  int kg = pt.get<int>("grid.kg", 2);
  assert(kg == 2);
  double tau_ini = pt.get<double>("adapt.tau_ini", 0.05);

  std::cout << "Compute a reference solution..." << std::endl;
  using namespace Dune::Functions::BasisFactory;
  auto feBasis = makeFlowBasis<2>(hostGrid->leafGridView());
  auto solution = run<Flow>(pt, tau_ini, feBasis, identityMap, "output2_ref.pvd", threadCount);
  auto positionBasis = Dune::Functions::subspaceBasis(feBasis, Dune::Indices::_0);
  auto X = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,3>>(positionBasis, solution);

  std::vector<double> errs, hs;
  for (int level = 0; level <= maxLevel; ++level) {
    std::cout << "Compute solution on level " << level << " ..." << std::endl;
    auto feBasis0 = makeFlowBasis<2>(hostGrid->levelGridView(level));
    auto solution0 = run<Flow>(pt, tau_ini,  feBasis0, identityMap, "output2_" + std::to_string(level) + ".pvd", threadCount);
    auto positionBasis0 = Dune::Functions::subspaceBasis(feBasis0, Dune::Indices::_0);
    auto X0 = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,3>>(positionBasis0, solution0);

    auto error0 = Dune::MeshDist::meanSquareError(X,X0,pt.sub("hausdorff"));
    std::cout << "error(" << level << ")  = " << error0  << std::endl;

    errs.push_back(error0);
    hs.push_back(Dune::BGN::gridSize(hostGrid->levelGridView(level)));
  }

  Dune::printErrorsClassic(std::cout, hs, {"dist(G,Gh)"}, {errs});
}

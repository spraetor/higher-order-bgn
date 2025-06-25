#pragma once

#include <cmath>
#include <iostream>

#include <dune/common/indices.hh>
#include <dune/common/parametertree.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/meshdist/hausdorffdistance.hh>

#include "flowbasis.hh"
#include "gridsize.hh"
#include "printerrors.hh"
#include "runner.hh"

namespace Dune::BGN {

// Compute the experimental order of convergence (EOC) for a flow problem
// where the error is defined in terms of the (symmetric) Hausdorff distance between
// the computed discrete surfaces
template <int kg, template<class> class Runner = ThreadRunner, class HostGrid, class InitialSurface, class Flow>
void eoc (const Dune::ParameterTree& pt, HostGrid& hostGrid, const InitialSurface& initialSurface, const Flow& flow, std::string outputBase = "output")
{
  using namespace Dune;

  int refinement_levels = pt.get<int>("grid.refinement_levels", 3);
  hostGrid.globalRefine(refinement_levels+2);

  double tau_ini = pt.get<double>("adapt.tau_ini", 0.05);

  std::cout << "Compute a reference solution..." << std::endl;
  auto feBasis = makeFlowBasis<kg>(hostGrid.leafGridView());
  double tau = tau_ini/std::pow(2, (kg+1)*(refinement_levels+2));
  auto runner = Runner<decltype(feBasis)>{feBasis, pr, tau};
  runner.init(initialSurface);
  runner.run(flow, outputBase + "_ref.pvd");
  auto positionBasis = Functions::subspaceBasis(feBasis, Indices::_0);
  auto X = Functions::makeDiscreteGlobalBasisFunction<FieldVector<double,3>>(positionBasis, runner.solution());

  std::vector<double> errs, hs;
  for (int level = 0; level <= hostGrid.maxLevel()-2; ++level) {
    std::cout << "Compute solution on level " << level << " ..." << std::endl;
    auto feBasis0 = makeFlowBasis<kg>(hostGrid.levelGridView(level));
    double tau_level = tau_ini/std::pow(2, (kg+1)*(level));
    auto runner0 = Runner<decltype(feBasis)>{feBasis0, pr, tau_level};
    runner0.init(initialSurface);
    runner0.run(flow, outputBase + "_" + std::to_string(level) + ".pvd");
    auto positionBasis0 = Functions::subspaceBasis(feBasis0, Indices::_0);
    auto X0 = Functions::makeDiscreteGlobalBasisFunction<FieldVector<double,3>>(positionBasis0, runner0.solution());

    auto error0 = MeshDist::meanSquareError(X,X0,pt.sub("hausdorff"));
    std::cout << "error(" << level << ")  = " << error0  << std::endl;

    errs.push_back(error0);
    hs.push_back(gridSize(hostGrid.levelGridView(level)));
  }

  printErrorsClassic(std::cout, hs, {"dist(G,Gh)"}, {errs});
}

} // end namepace Dune::BGN
#pragma once

#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

namespace Dune::BGN {

// Create a function-space basis [V^3 x V] suitable for the geometric flow problem
template <int kg = 2, class GridView>
auto makeFlowBasis(GridView const& gridView)
{
  using namespace Dune::Functions::BasisFactory;
  return makeBasis(gridView,
      composite(power<GridView::dimensionworld>(lagrange<kg>(), flatInterleaved()),
                lagrange<kg>(),
                flatLexicographic()));
}

} // end namespace Dune::BGN

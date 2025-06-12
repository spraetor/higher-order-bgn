#pragma once

#include <dune/curvedgeometry/localfunctiongeometry.hh>
#include <dune/geometry/quadraturerules.hh>

namespace Dune::BGN {

template <class GridView, class GridFct>
double surface(GridView const& gridView, GridFct const& X, int quadOrder = 10)
{
  auto X_e = localFunction(X);

  double vol = 0;
  for (auto const& e : elements(gridView)) {
    X_e.bind(e);
    auto geometry = Dune::LocalFunctionGeometry{referenceElement(e), X_e};
    const auto& quad = Dune::QuadratureRules<double, GridView::dimension>::rule(e.type(), quadOrder);
    vol += geometry.volume(quad);
  }
  return vol;
}

} // end namespace Dune::BGN

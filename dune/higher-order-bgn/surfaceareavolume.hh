#pragma once

#include <dune/curvedgeometry/localfunctiongeometry.hh>
#include <dune/geometry/quadraturerules.hh>

namespace Dune::BGN {

template <class GridView, class GridFct>
auto surfaceAreaVolume(GridView const& gridView, GridFct const& X, int quadOrder = 10)
{
  auto X_e = localFunction(X);

  double area = 0, vol = 0;
  for (auto const& e : elements(gridView)) {
    X_e.bind(e);
    auto geometry = Dune::LocalFunctionGeometry{referenceElement(e), X_e};
    const auto& quad = Dune::QuadratureRules<double, GridView::dimension>::rule(e.type(), quadOrder);
    for (auto const& [x,w] : quad) {
      const auto dS = geometry.integrationElement(x) * w;
      auto const& n = geometry.normal(x);

      // surface_area = int_{gridview} dx
      area += dS;

      auto fAtQP = X_e(x);
      for (std::size_t i = 0; i < fAtQP.size(); ++i) {
        // volume = int_{gridview} 1/3*X_e*normal dx
        vol += 1.0/fAtQP.size() * fAtQP[i] * n[i] * dS;
      }
    }
  }
  return std::tuple{area,vol};
}

} // end namespace Dune::BGN

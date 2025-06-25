#pragma once

#include <algorithm>

namespace Dune::BGN {

template <class GridView>
double gridSize(GridView const& gridView)
{
  using std::max;
  double h = 0;
  for (auto const& e : elements(gridView, Dune::Partitions::interior)) {
    for (unsigned int i = 0; i < e.subEntities(GridView::dimension-1); ++i)
      h = max(h, e.template subEntity<GridView::dimension - 1>(i).geometry().volume());
  }
  return h;
}

} // end namespace Dune::BGN

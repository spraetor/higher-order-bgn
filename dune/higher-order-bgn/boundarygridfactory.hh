#ifndef DUNE_BOUNDARY_GRIDFACTORY_HH
#define DUNE_BOUNDARY_GRIDFACTORY_HH

#include <cmath>
#include <limits>
#include <map>
#include <vector>

#include <dune/grid/common/gridfactory.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/utility/typefromvertexcount.hh>

namespace Dune
{
  template <class GridType>
  class BoundaryGridFactory
  {
    using ctype = typename GridType::ctype;

    enum { dim = GridType::dimension };
    enum { dimworld = GridType::dimensionworld };

    struct CoordLess
    {
      template <class T, int N>
      bool operator() (FieldVector<T,N> const& lhs, FieldVector<T,N> const& rhs) const
      {
        for (int i = 0; i < N; ++i) {
          using std::abs;
          if (abs(lhs[i] - rhs[i]) < std::numeric_limits<T>::epsilon())
            continue;
          return lhs[i] < rhs[i];
        }
        return false;
      }
    };

  public:
    /// Create a grid from the boundary intersections of a hostgrid that are accepted by the indicator
    template <class HostGridView, class Indicator>
    static void createBoundaryGrid (GridFactory<GridType>& factory,
                                    HostGridView const& hostGridView,
                                    Indicator indicator)
    {
      unsigned int idx = 0;
      std::map<FieldVector<ctype, dimworld>, unsigned int, CoordLess> vertices;

      for (auto const& e : elements(hostGridView)) {
        if (not e.hasBoundaryIntersections())
          continue;

        for (auto const& is : intersections(hostGridView, e)) {
          if (not is.boundary())
            continue;

          if (indicator(is)) {
            auto geo = is.geometry();
            std::vector<unsigned int> indices(geo.corners());
            for (int i = 0; i < geo.corners(); ++i) {
              auto it = vertices.emplace(geo.corner(i), idx);
              if (it.second) {
                ++idx;
                factory.insertVertex(geo.corner(i));
              }
              indices[i] = it.first->second;
            }
            factory.insertElement(geometryTypeFromVertexCount(dim, indices.size()), indices);
          }
        }
      }
    }

    /// Create a grid from all boundary intersections of a hostgrid
    template <class HostGridView>
    static void createBoundaryGrid (GridFactory<GridType>& factory,
                                    HostGridView const& hostGridView)
    {
      createBoundaryGrid(factory, hostGridView, [](auto const& is) { return true; });
    }

    /// Create a grid from the boundary intersections of a hostgrid that are accepted by the indicator
    template <class HostGridView, class Indicator>
    static std::unique_ptr<GridType> createBoundaryGrid (HostGridView const& hostGridView,
                                                         Indicator indicator)
    {
      GridFactory<GridType> factory;
      createBoundaryGrid(factory, hostGridView, indicator);
      return std::unique_ptr<GridType>(factory.createGrid());
    }

    /// Create a grid from all boundary intersections of a hostgrid
    template <class HostGridView>
    static std::unique_ptr<GridType> createBoundaryGrid (HostGridView const& hostGridView)
    {
      GridFactory<GridType> factory;
      createBoundaryGrid(factory, hostGridView, [](auto const& is) { return true; });
      return std::unique_ptr<GridType>(factory.createGrid());
    }
  };

} // end namespace Dune

#endif // DUNE_BOUNDARY_GRIDFACTORY_HH

#ifndef DUNE_VTK_DATACOLLECTORS_LOCALFUNCTIONGEOMETRYDATACOLLECTOR_HH
#define DUNE_VTK_DATACOLLECTORS_LOCALFUNCTIONGEOMETRYDATACOLLECTOR_HH

#include <vector>

#include <dune/common/referencehelper.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/vtk/types.hh>
#include <dune/vtk/datacollectors/unstructureddatacollector.hh>

namespace Dune::Vtk {

/// \brief DataCollector for dune-vtk writers extracting the parametrized geometry from
/// a local-to-global coordinate mapping and writing quadratic VTK elements
/**
 * To write curved elements, a \ref CurvedGeometry is constructed on-the-fly from the given
 * mappings in the constructor.
 *
 * \tparam GridView  The grid-view to write
 * \tparam Geometry  The CurvedGeometry to use for element parametrization
 **/
template <class GridView, class GF>
class LocalFunctionGeometryDataCollector
    : public UnstructuredDataCollectorInterface<GridView,
        LocalFunctionGeometryDataCollector<GridView,GF>, Partitions::All>
{
  using Self = LocalFunctionGeometryDataCollector;
  using Super = Vtk::UnstructuredDataCollectorInterface<GridView, Self, Partitions::All>;
  using GridFunction = GF;

  static auto nodeMapper (const GridView& gridView, Vtk::CellType::Parametrization param)
  {
    auto layout = [param](Dune::GeometryType gt, int dimgrid) -> unsigned int {
      auto cellType = Dune::Vtk::CellType{gt,param};
      return cellType.layout(gt.dim());
    };
    return Dune::MultipleCodimMultipleGeomTypeMapper(gridView, layout);
  }

  GridFunction gf_;
  Vtk::CellType::Parametrization parametrization_;
  Dune::MultipleCodimMultipleGeomTypeMapper<GridView> mapper_;

public:
  using Super::dim;
  using Super::partition; // NOTE: cell-type data-collector currently implemented for the All partition only
  using Super::gridView;

public:
  LocalFunctionGeometryDataCollector (GridView const& gridView, const GridFunction& gf, Vtk::CellType::Parametrization param = Vtk::CellType::LINEAR)
    : Super(gridView)
    , gf_(gf)
    , parametrization_(param)
    , mapper_(nodeMapper(Super::gridView(), param))
  {}

  /// Update the mapper
  void updateImpl ()
  {
    mapper_.update(Super::gridView());
  }

  /// Return number of vertices + number of edge
  std::uint64_t numPointsImpl () const
  {
    return mapper_.size();
  }

  /// Return a vector of point coordinates.
  /**
  * The vector of point coordinates is composed of vertex coordinates first and second
  * edge center coordinates.
  **/
  template <class T>
  std::vector<T> pointsImpl () const
  {
    std::vector<T> data(this->numPoints() * 3);
    auto lf = localFunction(Dune::resolveRef(gf_));
    for (auto const& c : elements(gridView(), partition)) {
      lf.bind(c);
      Vtk::CellType cellType(c.type(), parametrization_);
      auto refElem = referenceElement<T,dim>(c.type());

      for (int codim = dim, j = 0; codim >= 0 && j < cellType.size(); --codim) {
        for (int i = 0; i < refElem.size(codim) && j < cellType.size(); ++i,++j) {
          std::int64_t idx = 3 * mapper_.subIndex(c,i,codim);
          auto v = lf(refElem.position(i,codim));
          for (std::size_t j = 0; j < v.size(); ++j)
            data[idx + j] = T(v[j]);
          for (std::size_t j = v.size(); j < 3u; ++j)
            data[idx + j] = T(0);
        }
      }
    }
    return data;
  }

  /// Return number of grid cells
  std::uint64_t numCellsImpl () const
  {
    return gridView().size(0);
  }

  /// \brief Return cell types, offsets, and connectivity. \see Cells
  /**
  * The cell connectivity is composed of cell vertices first and second cell edges,
  * where the indices are grouped [vertex-indices..., (#vertices)+edge-indices...]
  **/
  Cells cellsImpl () const
  {
    Cells cells;
    cells.connectivity.reserve(this->numPoints());
    cells.offsets.reserve(this->numCells());
    cells.types.reserve(this->numCells());

    std::int64_t old_o = 0;
    for (auto const& c : elements(gridView(), partition)) {
      Vtk::CellType cellType(c.type(), parametrization_);
      auto refElem = referenceElement(c);

      for (int codim = dim, j = 0; codim >= 0 && j < cellType.size(); --codim) {
        for (int i = 0; i < refElem.size(codim) && j < cellType.size(); ++i,++j) {
          int k = cellType.permutation(j);
          cells.connectivity.push_back(mapper_.subIndex(c,k,codim));
        }
      }
      cells.offsets.push_back(old_o += cellType.size());
      cells.types.push_back(cellType.type());
    }
    return cells;
  }

  /// Evaluate the `fct` at element vertices and edge centers in the same order as the point coords.
  template <class T, class GlobalFunction>
  std::vector<T> pointDataImpl (GlobalFunction const& fct) const
  {
    std::vector<T> data(this->numPoints() * fct.numComponents());
    auto localFct = localFunction(fct);
    for (auto const& c : elements(gridView(), partition)) {
      localFct.bind(c);
      Vtk::CellType cellType{c.type(), parametrization_};
      auto refElem = referenceElement(c);

      for (int codim = dim, j = 0; codim >= 0 && j < cellType.size(); --codim) {
        for (int i = 0; i < refElem.size(codim) && j < cellType.size(); ++i,++j) {
          int k = cellType.permutation(j);
          std::int64_t idx = fct.numComponents() * mapper_.subIndex(c,k,codim);
          auto v = refElem.position(k,codim);
          for (int comp = 0; comp < fct.numComponents(); ++comp)
            data[idx + comp] = T(localFct.evaluate(comp, v));
        }
      }
      localFct.unbind();
    }
    return data;
  }
};

} // end namespace Dune::Vtk

#endif // DUNE_VTK_DATACOLLECTORS_LOCALFUNCTIONGEOMETRYDATACOLLECTOR_HH

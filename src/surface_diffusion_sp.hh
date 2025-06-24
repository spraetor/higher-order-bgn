#pragma once

#include <iostream>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/curvedgeometry/localfunctiongeometry.hh>
#include <dune/geometry/quadraturerules.hh>


namespace Dune::BGN {

template <class LocalView, class GridFct1, class GridFct2, class GridFct3>
class SDSPLocalAssembler
{
  using GridView = typename LocalView::GridView;
  using Element = typename LocalView::Element;
  using LocalFct1 = std::decay_t<decltype(localFunction(std::declval<Dune::ResolveRef_t<GridFct1> const&>()))>;
  using LocalFct2 = std::decay_t<decltype(localFunction(std::declval<Dune::ResolveRef_t<GridFct2> const&>()))>;
  using LocalFct3 = std::decay_t<decltype(localFunction(std::declval<Dune::ResolveRef_t<GridFct3> const&>()))>;

  static const int dim = Element::dimension;
  static const int dow = GridView::dimensionworld;


  static auto perp (FieldMatrix<double,2,3> const& J)
  {
    return FieldVector<double,3>{
      J[0][1] * J[1][2] - J[0][2] * J[1][1],
      J[0][2] * J[1][0] - J[0][0] * J[1][2],
      J[0][0] * J[1][1] - J[0][1] * J[1][0]};
   }

public:
  template <class Basis>
  SDSPLocalAssembler (const Basis&, const GridFct1& X, const GridFct2& Xmid, const GridFct2& Xiter, int quadOrder, double tau)
    : X_(X)
    , Xmid_(Xmid)
    , Xiter_(Xiter)
    // , X_e(localFunction(Dune::resolveRef(X_)))
    , quadOrder_(quadOrder)
    , tau_(tau)
  {}

  void bindLocalViews (const LocalView& testLocalView, const LocalView& trialLocalView)
  {
    testLocalView_ = &testLocalView;
    trialLocalView_ = &trialLocalView;
    X_e.emplace(localFunction(Dune::resolveRef(X_)));
    Xmid_e.emplace(localFunction(Dune::resolveRef(Xmid_)));
    Xiter_e.emplace(localFunction(Dune::resolveRef(Xiter_)));
  }

  void bindElement (const Element& element)
  {
    X_e->bind(element);
    Xmid_e->bind(element);
    Xiter_e->bind(element);
  }

  template <class LocalPattern>
  void assembleElementMatrixPattern (const Element& e, LocalPattern& p) noexcept
  {
    using namespace Dune::Indices;
    const auto& node = testLocalView_->tree();
    const auto& localFE = node.child(_1).finiteElement();

    for (auto i : Dune::range(localFE.size())) {
      for (auto j : Dune::range(localFE.size())) {
        auto row = node.child(_1).localIndex(i);
        auto col = node.child(_1).localIndex(j);
        p.add(row,col);

        for (auto k : Dune::range(node.child(_0).degree())) {
          auto row_k = node.child(_0).child(k).localIndex(i);
          auto col_k = node.child(_0).child(k).localIndex(j);
          p.add(row,col_k);
          p.add(row_k,col);
          p.add(row_k,col_k);
        }
      }
    }
  }

  template <class LocalMatrix>
  void assembleElementMatrix (const Element& e, LocalMatrix& localMatrix)
  {
    auto const& localFEmid = testLocalView_->tree().child(Indices::_0).child(0).finiteElement();

    auto geometry = Dune::LocalFunctionGeometry{referenceElement(e), *X_e};
    auto geometry_mid = Dune::ParametrizedGeometry{referenceElement(e), localFEmid, *Xmid_e};
    auto geometry_iter = Dune::LocalFunctionGeometry{referenceElement(e), *Xiter_e};

    using namespace Dune::Indices;

    const auto& node = testLocalView_->tree();
    const auto& localFE = node.child(_1).finiteElement();
    assert(node.child(_0).degree() == dow);

    values_.resize(localFE.size());
    jacobians_.resize(localFE.size());
    refJacobians_.resize(localFE.size());

    const auto& quad = Dune::QuadratureRules<double, dim>::rule(e.type(), quadOrder_);
    for (const auto& [x,w] : quad)
    {
      const auto integrationElement = geometry.integrationElement(x);
      const auto dx = integrationElement * w;
      const auto Jit = geometry.jacobianInverse(x);

      auto const J1 = geometry.jacobianTransposed(qp.position());
      auto const J2 = geometry_mid.jacobianTransposed(qp.position());
      auto const J3 = geometry_iter.jacobianTransposed(qp.position());

      const auto n_Picard = (perp(J1) + 4* perp(J2) + perp(J3)) / (6 * integrationElement);

      localFE.localBasis().evaluateFunction(x, values_);
      localFE.localBasis().evaluateJacobian(x, refJacobians_);
      for (auto i : Dune::range(localFE.size()))
        jacobians_[i] = refJacobians_[i] * Jit;

      for (auto i : Dune::range(localFE.size())) {
        for (auto j : Dune::range(localFE.size())) {
          auto row = node.child(_1).localIndex(i);
          auto col = node.child(_1).localIndex(j);

          // tau * (grad(v), grad(w))
          const auto laplace = jacobians_[i][0].dot(jacobians_[j][0]);
          localMatrix[row][col] -= tau_ * laplace * dx;

          for (auto k : Dune::range(node.child(_0).degree())) {
            auto row_k = node.child(_0).child(k).localIndex(i);
            auto col_k = node.child(_0).child(k).localIndex(j);

            // (v, normal*w)
            const auto mass_normal = values_[i] * values_[j] * n_Picard[k];
            localMatrix[row][col_k] += mass_normal * dx;
            localMatrix[row_k][col] += mass_normal * dx;

            // (grad(v), grad(w))
            localMatrix[row_k][col_k] += laplace * dx;
          }
        }
      }
    }
  }

  void bindLocalView (const LocalView& testLocalView)
  {
    testLocalView_ = &testLocalView;
    X_e.emplace(localFunction(Dune::resolveRef(X_)));
    Xmid_e.emplace(localFunction(Dune::resolveRef(Xmid_)));
    Xiter_e.emplace(localFunction(Dune::resolveRef(Xiter_)));
  }

  template <class LocalVector>
  void assembleElementVector (const Element& e, LocalVector& localVector)
  {
    using namespace Dune::Indices;
    auto geometry = Dune::LocalFunctionGeometry{referenceElement(e), *X_e};
    auto geometry_mid = Dune::LocalFunctionGeometry{referenceElement(e), *Xmid_e};
    auto geometry_iter = Dune::LocalFunctionGeometry{referenceElement(e), *Xiter_e};

    const auto& node = testLocalView_->tree();
    const auto& localFE = node.child(_1).finiteElement();
    values_.resize(localFE.size());

    const auto& quad = Dune::QuadratureRules<double, dim>::rule(e.type(), quadOrder_);
    for (const auto& [x,w] : quad)
    {
      const auto integrationElement = geometry.integrationElement(x);
      const auto dx = integrationElement * w;

      auto const J1 = geometry.jacobianTransposed(qp.position());
      auto const J2 = geometry_mid.jacobianTransposed(qp.position());
      auto const J3 = geometry_iter.jacobianTransposed(qp.position());

      const auto n_Picard = (perp(J1) + 4* perp(J2) + perp(J3)) / (6 * integrationElement);

      localFE.localBasis().evaluateFunction(x, values_);

      auto X_n = (*X_e)(x).dot(n_Picard);
      for (auto i : Dune::range(localFE.size())) {
        auto row = node.child(_1).localIndex(i);
        localVector[row] += values_[i] * X_n * dx;
      }
    }
  }

private:
  GridFct1 X_;
  GridFct2 Xmin_;
  GridFct3 Xiter_;
  std::optional<LocalFct1> X_e;
  std::optional<LocalFct2> Xmid_e;
  std::optional<LocalFct3> Xiter_e;
  int quadOrder_;
  double tau_;

  const LocalView* testLocalView_ = nullptr;
  const LocalView* trialLocalView_ = nullptr;
  std::vector<Dune::FieldMatrix<double,1,dim>> refJacobians_;
  std::vector<Dune::FieldMatrix<double,1,dow>> jacobians_;
  std::vector<Dune::FieldVector<double,1>> values_;
};

template <class Basis, class GridFct1, class GridFct2, class GridFct3>
SDSPLocalAssembler(const Basis&, const GridFct1&, const GridFct2&, const GridFct3&, int, double)
  -> SDSPLocalAssembler<typename Basis::LocalView, GridFct1, GridFct2, GridFct2>;

struct SurfaceDiffusionSP
{
  template <class Basis, class GridFct, class GridFct2>
  auto operator() (const Basis& basis, const GridFct1& X, const GRidFct3& Xiter, int quadOrder, double tau) const
  {
    auto Xmid = Dune::Functions::makeComposedGridFunction([](auto X_x,auto X_iter_x){
      return (X_x + X_iter_x)/2.0;}, X, Xiter);
    return SDSPLocalAssembler{basis, X, Xmid, Xiter, quadOrder,tau};
  }
};

} // end namespace Dune::BGN

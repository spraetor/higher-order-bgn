#pragma once

#include <iostream>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/curvedgeometry/localfunctiongeometry.hh>
#include <dune/geometry/quadraturerules.hh>


namespace Dune::BGN {

template <class LocalView, class GridFct>
class SDLocalAssembler
{
  using GridView = typename LocalView::GridView;
  using Element = typename LocalView::Element;
  using LocalFct = std::decay_t<decltype(localFunction(std::declval<Dune::ResolveRef_t<GridFct> const&>()))>;

  static const int dim = Element::dimension;
  static const int dow = GridView::dimensionworld;

public:
  template <class Basis>
  SDLocalAssembler (const Basis&, const GridFct& X, int quadOrder, double tau)
    : X_(X)
    // , X_e(localFunction(Dune::resolveRef(X_)))
    , quadOrder_(quadOrder)
    , tau_(tau)
  {}

  void bindLocalViews (const LocalView& testLocalView, const LocalView& trialLocalView)
  {
    testLocalView_ = &testLocalView;
    trialLocalView_ = &trialLocalView;
    X_e.emplace(localFunction(Dune::resolveRef(X_)));
  }

  void bindElement (const Element& element)
  {
    X_e->bind(element);
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
    auto geometry = Dune::LocalFunctionGeometry{referenceElement(e), *X_e};

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
      const auto dx = geometry.integrationElement(x) * w;
      const auto Jit = geometry.jacobianInverse(x);
      const auto n = geometry.normal(x);

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
            const auto mass_normal = values_[i] * values_[j] * n[k];
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
  }

  template <class LocalVector>
  void assembleElementVector (const Element& e, LocalVector& localVector)
  {
    using namespace Dune::Indices;
    auto geometry = Dune::LocalFunctionGeometry{referenceElement(e), *X_e};

    const auto& node = testLocalView_->tree();
    const auto& localFE = node.child(_1).finiteElement();
    values_.resize(localFE.size());

    const auto& quad = Dune::QuadratureRules<double, dim>::rule(e.type(), quadOrder_);
    for (const auto& [x,w] : quad)
    {
      const auto dx = geometry.integrationElement(x) * w;
      const auto n = geometry.normal(x);

      localFE.localBasis().evaluateFunction(x, values_);

      auto X_n = (*X_e)(x).dot(n);
      for (auto i : Dune::range(localFE.size())) {
        auto row = node.child(_1).localIndex(i);
        localVector[row] += values_[i] * X_n * dx;
      }
    }
  }

private:
  GridFct X_;
  std::optional<LocalFct> X_e;
  int quadOrder_;
  double tau_;

  const LocalView* testLocalView_ = nullptr;
  const LocalView* trialLocalView_ = nullptr;
  std::vector<Dune::FieldMatrix<double,1,dim>> refJacobians_;
  std::vector<Dune::FieldMatrix<double,1,dow>> jacobians_;
  std::vector<Dune::FieldVector<double,1>> values_;
};

template <class Basis, class GridFct>
SDLocalAssembler(const Basis&, const GridFct& Xh, int, double)
  -> SDLocalAssembler<typename Basis::LocalView, GridFct>;

struct SurfaceDiffusion
{
  template <class Basis, class GridFct>
  auto operator() (const Basis& basis, const GridFct& X, int quadOrder, double tau) const
  {
    return SDLocalAssembler{basis,X,quadOrder,tau};
  }
};

} // end namespace Dune::BGN

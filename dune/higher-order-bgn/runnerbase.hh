#pragma once

#include <cmath>
#include <iostream>
#include <thread>

#include <dune/common/exceptions.hh>
#include <dune/common/indices.hh>
#include <dune/common/parametertree.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>

namespace Dune::BGN {

template <class B, class A>
struct RunnerBase : public A
{
  using Basis = B;
  using Assembler = A;

  using GridView = typename Basis::GridView;
  static constexpr int dow = GridView::dimensionworld;

  using Matrix = Dune::BCRSMatrix<double>;
  using Vector = Dune::BlockVector<double>;

  RunnerBase (const Basis& feBasis, const Dune::ParameterTree& pt, double tau)
    : Assembler(feBasis, pt.sub("assembler"))
    , feBasis_(feBasis)
    , pt_(pt)
    , tau_(tau)
  {}

  RunnerBase (const Basis& feBasis, const Dune::ParameterTree& pt)
    : RunnerBase{feBasis, pt, pt.get<double>("adapt.timestep")}
  {}

  RunnerBase (const RunnerBase& other)
    : RunnerBase{other.feBasis_, other.pt_, other.tau_}
  {}

  // compute an initial solution
  template <class InitialSurface>
  void init (InitialSurface const& initialSurface)
  {
    auto positionBasis = Dune::Functions::subspaceBasis(feBasis_, Dune::Indices::_0);

    // interpolate the initial surface parametrization into the solution vector
    Dune::Functions::interpolate(positionBasis, solution_, initialSurface);
  }

  template <class Flow>
  void run (const Flow& flow, std::string outputFileName = "")
  {
    DUNE_THROW(Dune::NotImplemented, "The run() method in this runner is not implemented");
  }

  Vector const& solution () const
  {
    return solution_;
  }

  void setTimestep (double tau)
  {
    tau_ = tau;
  }

  double timestep () const
  {
    return tau_;
  }

protected:
  const Basis& feBasis_;
  const Dune::ParameterTree& pt_;
  Vector solution_;
  double tau_;
};

} // end namespace Dune::BGN

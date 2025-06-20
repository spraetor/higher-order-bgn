#pragma once

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/curvedgrid/geometries/implicitsurface.hh>

namespace Dune::BGN {

template <class T = double>
class EllipseProjection
{
  using Domain = FieldVector<T,2>;
  using Jacobian = FieldMatrix<T,2,2>;


  // Implicit function representation of the ellipsoid surface
  struct Phi
  {
    T a_, b_;

    // phi(x,y) = (x/a)^2 + (y/b)^2 = 1
    T operator() (const FieldVector<T,2>& x) const
    {
      return x[0]*x[0]/(a_*a_) + x[1]*x[1]/(b_*b_) - 1;
    }

    // grad(phi)
    friend auto derivative (Phi phi)
    {
      return [a=phi.a_,b=phi.b_](const Domain& x) -> Domain
      {
        return { T(2*x[0]/(a*a)), T(2*x[1]/(b*b)) };
      };
    }
  };

public:
  /// \brief Constructor of ellipse by major axes
  EllipseProjection (T a, T b)
    : a_(a)
    , b_(b)
    , implicit_(Phi{a,b},100)
  {}

  /// \brief project the coordinate to the ellipsoid
  Domain operator() (const Domain& x) const
  {
    return implicit_(x);
    }

  private:
    T a_, b_;
    SimpleImplicitSurfaceProjection<Phi> implicit_;
  };


} // end namespace Dune::BGN

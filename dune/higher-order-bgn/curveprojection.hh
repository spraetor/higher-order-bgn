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
  Dune::SimpleImplicitSurfaceProjection<Phi> implicit_;
};


template <class T = double>
class Circle
{
  using Domain = FieldVector<T,2>;
  using Jacobian = FieldMatrix<T,2,2>;

  // Implicit function representation of the ellipsoid surface
  struct Phi
  {
    T r_;

    // phi(x,y) = (x/r)^2 + (y/r)^2 = 1
    T operator() (const FieldVector<T,2>& x) const
    {
      return x[0]*x[0]/(r_*r_) + x[1]*x[1]/(r_*r_) - 1;
    }

    // grad(phi)
    friend auto derivative (Phi phi)
    {
      return [r=phi.r_](const Domain& x) -> Domain
      {
        return { T(2*x[0]/(r*r)), T(2*x[1]/(r*r)) };
      };
    }
  };

public:
  /// \brief Constructor of ellipse by major axes
  Circle (T r)
    : r_(r)
    , implicit_(Phi{r},100)
  {}

  /// \brief project the coordinate to the ellipsoid
  Domain operator() (const Domain& x) const
  {
    return implicit_(x);
  }

private:
  T r_;
  Dune::SimpleImplicitSurfaceProjection<Phi> implicit_;
};


class PerturbedCircle
{
public:
  PerturbedCircle(double radius, double twist, double period)
    : radius_(radius)
    , twist_(twist)
    , period_(period)
  {}

  FieldVector<double,2> operator()(FieldVector<double,2> const& x) const
  {
    using std::sqrt;
    using std::atan2;
    double alpha = atan2(x[1], x[0]);
    return x * (radius(alpha)/x.two_norm());
  }

private:
  double radius(double alpha) const
  {
    using std::cos;
    using std::sin;
    double factor = twist_ +  (cos(period_*alpha));
    return radius_ * factor;
  }

private:
  double radius_;
  double twist_;
  double period_;
};

} // end namespace Dune::BGN

#pragma once

#include <cmath>
#include <iostream>
#include <thread>
#include <vector>

#include <dune/assembler/backends/istlvectorbackend.hh>
#include <dune/assembler/backends/istlmatrixbackend.hh>
#include <dune/common/indices.hh>
#include <dune/common/parametertree.hh>
#include <dune/curvedgeometry/localfunctiongeometry.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/umfpack.hh>
#include <dune/vtk/pvdwriter.hh>
#include <dune/vtk/vtkwriter.hh>

#include "localfunctiongeometrydatacollector.hh"
#include "runner.hh"
#include "surfaceareavolume.hh"
#include "umfpack2.hh"


namespace Dune::BGN {

template <class B, class A = ThreadAssembler<B>>
struct PicardRunner : public RunnerBase<B,A>
{
  using Basis = B;
  using Assembler = A;
  using Super = RunnerBase<B,A>;

  using GridView = typename Basis::GridView;
  static constexpr int dow = GridView::dimensionworld;

  using Matrix = Dune::BCRSMatrix<double>;
  using Vector = Dune::BlockVector<double>;

  PicardRunner (const Basis& feBasis, const Dune::ParameterTree& pt, double tau)
    : Super{feBasis, pt, tau}
  {}

  PicardRunner (const Basis& feBasis, const Dune::ParameterTree& pt)
    : Super{feBasis, pt}
  {}

  template <class Flow>
  void run (const Flow& flow, std::string outputFileName = "")
  {
    using namespace Dune;
    Timer timer;

    Vector solution_iter;

    auto positionBasis = Dune::Functions::subspaceBasis(feBasis_, Indices::_0);
    auto Xh = Dune::Functions::makeDiscreteGlobalBasisFunction<FieldVector<double,dow>>(positionBasis, solution_);
    auto Xh_iter = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,dow>>(positionBasis, solution_iter);

    // define a (parallel) assembler
    auto matrix = Matrix();
    auto matrixBackend = Dune::Assembler::ISTLMatrixBackend(matrix);
    auto rhs = Vector();
    auto rhsBackend = Dune::Assembler::ISTLVectorBackend(rhs);

    int quadOrder = pt_.template get<double>("solution.quad_order",10);
    auto localAssembler = flow(feBasis_, std::cref(Xh), std::cref(Xh_iter), quadOrder, tau_);

    auto patternBuilder = matrixBackend.patternBuilder();
    patternBuilder.resize(feBasis_, feBasis_);
    this->assembleMatrixPattern(localAssembler, patternBuilder);
    patternBuilder.setupMatrix();

    // create linear solver
    auto solver = BGN::UMFPack<Matrix>{};
    solver.setVerbosity(pt_.template get<int>("solver.verbose"));

    rhsBackend.resize(feBasis_);

    double startTime = pt_.template get<double>("adapt.start_time", 0.0);
    double endTime = pt_.template get<double>("adapt.end_time", 0.03);
    if (outputFileName.empty())
      outputFileName = pt_.template get<std::string>("output.filename", "output.pvd");

    using std::sqrt;
    double tol = sqrt(std::numeric_limits<double>::epsilon());
    double iter_tol = (pt_.template get<double>("adapt.iter_tol",1e-12));
    int maxIter = (pt_.template get<double>("adapt.maxIter",100));

    Vtk::LocalFunctionGeometryDataCollector dataCollector{feBasis_.gridView(), Xh, Vtk::CellType::QUADRATIC};
    Vtk::UnstructuredGridWriter writer{dataCollector};
    Vtk::PvdWriter pvdWriter{writer};

    writer.write(outputFileName);
    pvdWriter.writeTimestep(startTime, outputFileName, "_piecefiles");

    std::vector<double> iteration_error_vector;
    std::vector<double> iteration_number_vector;

    solution_iter = solution_;

    double t = startTime;
    int nSteps = int(std::ceil((endTime - startTime)/tau_));
    for (int step = 0; t + tau_ < endTime + tol; ++step, t+= tau_)
    {
      std::cout << step << "/" << nSteps;
      auto [area,volume] = BGN::surfaceAreaVolume(feBasis_.gridView(), Xh);
      std::cout << " area = " << area;
      std::cout << " volume = " << volume;
      std::cout << std::endl;
      if ((step+1) % std::max(1,nSteps/100) == 0) {
        pvdWriter.writeTimestep(t, outputFileName, "_piecefiles");
      }

      double iter_error = 1e10;

      int iter_number = 0;
      for (; iter_number < maxIter && iter_error > iter_tol; ++iter_number)
      {
        auto solution_iter_last = solution_iter;

        // assemble matrix and vector
        matrixBackend.setZero();
        this->assembleMatrixEntries(localAssembler, matrixBackend);
        rhsBackend.setZero();
        this->assembleVectorEntries(localAssembler, rhsBackend);

        // solve the linear system
        InverseOperatorResult statistics;
        solver.setMatrix(matrix, {}, step>0);
        solver.apply(solution_iter, rhs, statistics);

        // compute iteration error
        solution_iter_last -= solution_iter;
        iter_error = solution_iter_last.two_norm();
        iteration_error_vector.push_back(iter_error);
      }
      iteration_number_vector.push_back(iter_number);
      solution_ = solution_iter;
    }
    pvdWriter.writeTimestep(endTime, outputFileName, "_piecefiles");

    std::cout << "elapsed time: " << timer.elapsed() << std::endl;
  }

protected:
  using Super::feBasis_;
  using Super::pt_;
  using Super::solution_;
  using Super::tau_;
};

} // end namespace Dune::BGN

#pragma once

#include <cmath>
#include <iostream>
#include <thread>

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

#include "assembler.hh"
#include "flowbasis.hh"
#include "localfunctiongeometrydatacollector.hh"
#include "runnerbase.hh"
#include "surfaceareavolume.hh"
#include "umfpack2.hh"


namespace Dune::BGN {

template <class B, class A>
struct Runner : public RunnerBase<B,A>
{
  using Basis = B;
  using Assembler = A;
  using Super = RunnerBase<B,A>;

  using GridView = typename Basis::GridView;
  static constexpr int dow = GridView::dimensionworld;

  using Matrix = Dune::BCRSMatrix<double>;
  using Vector = Dune::BlockVector<double>;

  Runner (const Basis& feBasis, const Dune::ParameterTree& pt, double tau)
    : RunnerBase<B,A>(feBasis, pt, tau)
  {}

  Runner (const Basis& feBasis, const Dune::ParameterTree& pt)
    : Runner{feBasis, pt, pt.get<double>("adapt.timestep")}
  {}

  template <class Flow>
  void run (const Flow& flow, std::string outputFileName = "")
  {
    using namespace Dune;
    Timer timer;

    auto positionBasis = Functions::subspaceBasis(feBasis_, Indices::_0);
    auto Xh = Functions::makeDiscreteGlobalBasisFunction<FieldVector<double,dow>>(positionBasis, solution_);

    // define a (parallel) assembler
    auto matrix = Matrix();
    auto matrixBackend = Dune::Assembler::ISTLMatrixBackend(matrix);
    auto rhs = Vector();
    auto rhsBackend = Dune::Assembler::ISTLVectorBackend(rhs);

    int quadOrder = pt_.template get<double>("solution.quad_order",10);
    auto localAssembler = flow(feBasis_, std::cref(Xh), quadOrder, tau_);

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

    Vtk::LocalFunctionGeometryDataCollector dataCollector{feBasis_.gridView(), Xh, Vtk::CellType::QUADRATIC};
    Vtk::UnstructuredGridWriter writer{dataCollector};
    Vtk::PvdWriter pvdWriter{writer};

    writer.write(outputFileName);
    pvdWriter.writeTimestep(startTime, outputFileName, "_piecefiles");

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

      // assemble matrix and vector
      matrixBackend.setZero();
      this->assembleMatrixEntries(localAssembler, matrixBackend);
      rhsBackend.setZero();
      this->assembleVectorEntries(localAssembler, rhsBackend);

      // solve the linear system
      InverseOperatorResult statistics;
      solver.setMatrix(matrix, {}, step>0);
      solver.apply(solution_, rhs, statistics);
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

template <class B>
using ThreadRunner = Runner<B, ThreadAssembler<B>>;

template <class B>
using DefaultRunner = Runner<B, DefaultAssembler<B>>;

template <int kg, template<class> class Runner = ThreadRunner, class HostGrid, class InitialSurface, class Flow>
void run (const Dune::ParameterTree& pt, HostGrid& hostGrid, const InitialSurface& initialSurface, const Flow& flow, std::string outputBase = "output")
{
  int refinement_levels = pt.get<int>("grid.refinement_levels", 3);
  hostGrid.globalRefine(refinement_levels);

  double tau = pt.get<double>("adapt.timestep", 0.05);

  std::cout << "Compute solution..." << std::endl;
  auto feBasis = makeFlowBasis<kg>(hostGrid.leafGridView());
  auto runner = Runner<decltype(feBasis)>{feBasis, pt, tau};
  runner.init(initialSurface);
  runner.run(flow, outputBase + ".pvd");
}

} // end namespace Dune::BGN

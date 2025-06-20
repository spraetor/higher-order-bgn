#pragma once

#include <cmath>
#include <iostream>
#include <thread>

#include <dune/assembler/backends/istlvectorbackend.hh>
#include <dune/assembler/backends/istlmatrixbackend.hh>
#include <dune/assembler/defaultglobalassembler.hh>
#include <dune/assembler/parallel/gridcoloring.hh>
#include <dune/assembler/parallel/coloredrangeexecutor.hh>
#include <dune/common/indices.hh>
#include <dune/common/parametertree.hh>
#include <dune/curvedgeometry/localfunctiongeometry.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/umfpack.hh>
#include <dune/meshdist/hausdorffdistance.hh>
#include <dune/vtk/pvdwriter.hh>
#include <dune/vtk/vtkwriter.hh>

#include "area.hh"
#include "gridsize.hh"
#include "localfunctiongeometrydatacollector.hh"
#include "printerrors.hh"
#include "umfpack2.hh"

// Create a function-space basis [V^3 x V] suitable for the geometric flow problem
template <int kg = 2, class GridView>
auto makeFlowBasis(GridView const& gridView)
{
  using namespace Dune::Functions::BasisFactory;
  return makeBasis(gridView,
      composite(power<GridView::dimensionworld>(lagrange<kg>(), flatInterleaved()),
                lagrange<kg>(),
                flatLexicographic()));
}


template <class Basis, class InitialSurface, class Flow>
auto run(const Dune::ParameterTree& pt, const Basis& feBasis, InitialSurface const& initialSurface, const Flow& flow, double tau, std::string outputFileName = "")
{
  using namespace Dune;

  Timer timer;
  using Matrix = BCRSMatrix<double>;
  using Vector = BlockVector<double>;

  using GridView = typename Basis::GridView;
  constexpr int dow = GridView::dimensionworld;

  // define a solution vector and associated discrete functions
  Vector solution;

  auto positionBasis = Functions::subspaceBasis(feBasis, Indices::_0);
  auto Xh = Functions::makeDiscreteGlobalBasisFunction<FieldVector<double,dow>>(positionBasis, solution);

  auto curvatureBasis = Functions::subspaceBasis(feBasis, Indices::_1);
  auto Hh = Functions::makeDiscreteGlobalBasisFunction<double>(curvatureBasis, solution);

  // interpolate the initial surface parametrization into the solution vector
  Functions::interpolate(positionBasis, solution, initialSurface);

  // number of threads to use
  std::size_t threadCount = pt.get<std::size_t>("solution.thread_count", std::thread::hardware_concurrency());

  // define a (parallel) assembler
  auto gridView = feBasis.gridView();
  const auto coloredRange = Assembler::Experimental::coloredElementRange(gridView);
  auto executor = Assembler::Experimental::ColoredRangeExecutor(coloredRange, threadCount);
  auto matrix = Matrix();
  auto matrixBackend = Assembler::ISTLMatrixBackend(matrix);
  auto rhs = Vector();
  auto rhsBackend = Assembler::ISTLVectorBackend(rhs);
  auto assembler = Assembler::Assembler(feBasis, feBasis, executor);

  int quadOrder = pt.get<double>("solution.quad_order",10);
  auto localAssembler = flow(feBasis, std::cref(Xh), quadOrder, tau);

  auto patternBuilder = matrixBackend.patternBuilder();
  patternBuilder.resize(feBasis, feBasis);
  assembler.assembleMatrixPattern(localAssembler, patternBuilder);
  patternBuilder.setupMatrix();

  // create linear solver
  auto solver = BGN::UMFPack<Matrix>{};
  solver.setVerbosity(pt.get<int>("solver.verbose"));

  rhsBackend.resize(feBasis);

  double startTime = pt.get<double>("adapt.start_time", 0.0);
  double endTime = pt.get<double>("adapt.end_time", 0.03);
  if (outputFileName.empty())
    outputFileName = pt.get<std::string>("output.filename", "output.pvd");

  using std::sqrt;
  double tol = sqrt(std::numeric_limits<double>::epsilon());

  Vtk::LocalFunctionGeometryDataCollector dataCollector{feBasis.gridView(), Xh, Vtk::CellType::QUADRATIC};
  Vtk::UnstructuredGridWriter writer{dataCollector};
  Vtk::PvdWriter pvdWriter{writer};

  writer.write(outputFileName);
  pvdWriter.writeTimestep(startTime, outputFileName, "_piecefiles");

  double t = startTime;

  int nSteps = int(std::ceil((endTime - startTime)/tau));
  for (int step = 0; t + tau < endTime + tol; ++step, t+= tau)
  {
    std::cout << step << "/" << nSteps;
    std::cout << " surface = " << BGN::surface(gridView, Xh);
    std::cout << std::endl;
    if ((step+1) % std::max(1,nSteps/100) == 0) {
      pvdWriter.writeTimestep(t, outputFileName, "_piecefiles");
    }

    // assemble matrix and vector
    matrixBackend.setZero();
    assembler.assembleMatrixEntries(localAssembler, matrixBackend);
    rhsBackend.setZero();
    assembler.assembleVectorEntries(localAssembler, rhsBackend);

    // solve the linear system
    InverseOperatorResult statistics;
    solver.setMatrix(matrix, {}, step>0);
    solver.apply(solution, rhs, statistics);
  }
  pvdWriter.writeTimestep(endTime, outputFileName, "_piecefiles");

  std::cout << "elapsed time: " << timer.elapsed() << std::endl;
  return solution;
}


template <int kg, class HostGrid, class InitialSurface, class Flow>
void run (const Dune::ParameterTree& pt, HostGrid& hostGrid, const InitialSurface& initialSurface, const Flow& flow, std::string outputBase = "output")
{
  using namespace Dune;

  int refinement_levels = pt.get<int>("grid.refinement_levels", 3);
  hostGrid.globalRefine(refinement_levels+2);

  double tau_ini = pt.get<double>("adapt.tau_ini", 0.05);

  std::cout << "Compute a reference solution..." << std::endl;
  using namespace Functions::BasisFactory;
  auto feBasis = makeFlowBasis<kg>(hostGrid.leafGridView());
  double tau = tau_ini/std::pow(2, (kg+1)*(refinement_levels+2));
  auto solution = run(pt, feBasis, initialSurface, flow, tau, outputBase + "_ref.pvd");
  auto positionBasis = Functions::subspaceBasis(feBasis, Indices::_0);
  auto X = Functions::makeDiscreteGlobalBasisFunction<FieldVector<double,3>>(positionBasis, solution);

  std::vector<double> errs, hs;
  for (int level = 0; level <= hostGrid.maxLevel()-2; ++level) {
    std::cout << "Compute solution on level " << level << " ..." << std::endl;
    auto feBasis0 = makeFlowBasis<kg>(hostGrid.levelGridView(level));
    double tau_level = tau_ini/std::pow(2, (kg+1)*(level));
    auto solution0 = run(pt, feBasis0, initialSurface, flow, tau_level, outputBase + "_" + std::to_string(level) + ".pvd");
    auto positionBasis0 = Functions::subspaceBasis(feBasis0, Indices::_0);
    auto X0 = Functions::makeDiscreteGlobalBasisFunction<FieldVector<double,3>>(positionBasis0, solution0);

    auto error0 = MeshDist::meanSquareError(X,X0,pt.sub("hausdorff"));
    std::cout << "error(" << level << ")  = " << error0  << std::endl;

    errs.push_back(error0);
    hs.push_back(BGN::gridSize(hostGrid.levelGridView(level)));
  }

  printErrorsClassic(std::cout, hs, {"dist(G,Gh)"}, {errs});
}
#pragma once

#include <dune/assembler/backends/istlvectorbackend.hh>
#include <dune/assembler/backends/istlmatrixbackend.hh>
#include <dune/assembler/defaultglobalassembler.hh>
#include <dune/assembler/parallel/gridcoloring.hh>
#include <dune/assembler/parallel/coloredrangeexecutor.hh>
#include <dune/common/parametertree.hh>
#include <dune/curvedgeometry/localfunctiongeometry.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/umfpack.hh>
#include <dune/vtk/pvdwriter.hh>
#include <dune/vtk/vtkwriter.hh>

#include "area.hh"
#include "localfunctiongeometrydatacollector.hh"
#include "umfpack2.hh"

template <class Flow, class Basis, class InitialSurface>
auto run(Dune::ParameterTree const& pt, const double& tau,  const Basis& feBasis, InitialSurface const& initialSurface, std::string outputFileName = "", std::size_t threadCount = 1)
{
  Dune::Timer timer;
  using Matrix = Dune::BCRSMatrix<double>;
  using Vector = Dune::BlockVector<double>;
  Vector solution;

  auto positionBasis = Dune::Functions::subspaceBasis(feBasis, Dune::Indices::_0);
  auto curvatureBasis = Dune::Functions::subspaceBasis(feBasis, Dune::Indices::_1);

  using GridView = typename Basis::GridView;
  auto Xh = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,GridView::dimensionworld>>(positionBasis, solution);
  auto Hh = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(curvatureBasis, solution);
  Dune::Functions::interpolate(positionBasis, solution, initialSurface);

  auto gridView = feBasis.gridView();
  const auto coloredRange = Dune::Assembler::Experimental::coloredElementRange(gridView);
  auto executor = Dune::Assembler::Experimental::ColoredRangeExecutor(coloredRange, threadCount);
  auto matrix = Matrix();
  auto matrixBackend = Dune::Assembler::ISTLMatrixBackend(matrix);
  auto rhs = Vector();
  auto rhsBackend = Dune::Assembler::ISTLVectorBackend(rhs);
  auto assembler = Dune::Assembler::Assembler(feBasis, feBasis, executor);

  int quadOrder = pt.get<double>("solution.quad_order",10);
  auto localAssembler = Flow{}(feBasis, std::cref(Xh), quadOrder, tau);

  auto patternBuilder = matrixBackend.patternBuilder();
  patternBuilder.resize(feBasis, feBasis);
  assembler.assembleMatrixPattern(localAssembler, patternBuilder);
  patternBuilder.setupMatrix();

  // create linear solver
  auto solver = Dune::BGN::UMFPack<Matrix>{};
  solver.setVerbosity(pt.get<int>("solver.verbose"));

  rhsBackend.resize(feBasis);

  double startTime = pt.get<double>("adapt.start_time", 0.0);
  double endTime = pt.get<double>("adapt.end_time", 0.03);
  if (outputFileName.empty())
    outputFileName = pt.get<std::string>("output.filename", "output.pvd");

  using std::sqrt;
  double tol = sqrt(std::numeric_limits<double>::epsilon());

  Dune::Vtk::LocalFunctionGeometryDataCollector dataCollector{feBasis.gridView(), Xh, Dune::Vtk::CellType::QUADRATIC};
  Dune::Vtk::UnstructuredGridWriter writer{dataCollector};
  Dune::Vtk::PvdWriter pvdWriter{writer};

  writer.write(outputFileName);
  pvdWriter.writeTimestep(startTime, outputFileName, "_piecefiles");

  double t = startTime;

  int nSteps = int(std::ceil((endTime - startTime)/tau));
  for (int step = 0; t + tau < endTime + tol; ++step, t+= tau)
  {
    std::cout << step << "/" << nSteps;
    std::cout << " surface = " << surface(gridView, Xh);
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
    Dune::InverseOperatorResult statistics;
    solver.setMatrix(matrix, {}, step>0);
    solver.apply(solution, rhs, statistics);
  }
  pvdWriter.writeTimestep(endTime, outputFileName, "_piecefiles");

  std::cout << "elapsed time: " << timer.elapsed() << std::endl;
  return solution;
}
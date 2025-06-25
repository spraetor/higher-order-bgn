#pragma once

#include <cmath>
#include <iostream>
#include <thread>

#include <dune/assembler/backends/istlvectorbackend.hh>
#include <dune/assembler/backends/istlmatrixbackend.hh>
#include <dune/assembler/defaultglobalassembler.hh>
#include <dune/assembler/parallel/gridcoloring.hh>
#include <dune/assembler/parallel/coloredrangeexecutor.hh>

namespace Dune::BGN {

template <class Basis>
struct DefaultAssembler : public Dune::Assembler::Assembler<Basis,Basis>
{
  DefaultAssembler (const Basis& feBasis, const Dune::ParameterStree& /*pt*/)
    : Dune::Assembler::Assembler<Basis,Basis>(feBasis,feBasis)
  {}
};


template <class Basis>
struct ThreadAssembler
{
  using GridView = typename Basis::GridView;
  using ColoredRange = decltype(Dune::Assembler::Experimental::coloredElementRange(std::declval<GridView>()));
  using Executor = Dune::Assembler::Experimental::ColoredRangeExecutor<ColoredRange>;
  using Assembler = Dune::Assembler::Assembler<Basis,Basis,Executor>;

public:
  ThreadAssembler (const Basis& feBasis, const Dune::ParameterStree& pt)
    : coloredRange_(Dune::Assembler::Experimental::coloredElementRange(feBasis_.gridView()))
    , executor_(Dune::Assembler::Experimental::ColoredRangeExecutor(coloredRange_,
        pt.get<std::size_t>("thread_count",std::thread::hardware_concurrency())))
    , assembler_(Dune::Assembler::Assembler(feBasis_, feBasis_, executor_))
  {}

  template <class LA, class PB>
  void assembleMatrixPattern (LA const& localAssembler, PB& patternBuilder)
  {
    assembler_.assembleMatrixPattern(localAssembler, patternBuilder);
  }

  template <class LA, class MB>
  void assembleMatrixEntries (LA const& localAssembler, MB& matrixBackend)
  {
    assembler_.assembleMatrixEntries(localAssembler, matrixBackend);
  }

  template <class LA, class VB>
  void assembleMatrixEntries (LA const& localAssembler, VB& rhsBackend)
  {
    assembler_.assembleVectorEntries(localAssembler, rhsBackend);
  }

private:
  using ColoredRange = decltype(Dune::Assembler::Experimental::coloredElementRange(std::declval<GridView>()));
  ColoredRange coloredRange_;

  using Executor = Dune::Assembler::Experimental::ColoredRangeExecutor<ColoredRange>;
  Executor executor_;

  using Assembler = Dune::Assembler::Assembler<Basis,Basis,Executor>;
  Assembler assembler_;
};

} // end namespace Dune::BGN

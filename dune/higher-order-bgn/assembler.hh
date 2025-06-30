#pragma once

#include <cmath>
#include <iostream>
#include <thread>

#include <dune/assembler/backends/istlvectorbackend.hh>
#include <dune/assembler/backends/istlmatrixbackend.hh>
#include <dune/assembler/defaultglobalassembler.hh>
#include <dune/assembler/parallel/gridcoloring.hh>
#include <dune/assembler/parallel/coloredrangeexecutor.hh>
#include <dune/common/parametertree.hh>

namespace Dune::BGN {

template <class Basis>
struct DefaultAssembler : public Dune::Assembler::Assembler<Basis,Basis>
{
  DefaultAssembler (const Basis& feBasis, const Dune::ParameterTree& /*pt*/)
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
  ThreadAssembler (const Basis& feBasis, const Dune::ParameterTree& pt)
    : coloredRange_(Dune::Assembler::Experimental::coloredElementRange(feBasis.gridView()))
    , executor_(Dune::Assembler::Experimental::ColoredRangeExecutor(coloredRange_,
        pt.get<std::size_t>("thread_count",std::thread::hardware_concurrency())))
    , assembler_(Dune::Assembler::Assembler(feBasis, feBasis, executor_))
  {}

  template <class LA, class PB>
  void assembleMatrixPattern (LA& localAssembler, PB& patternBuilder)
  {
    assembler_.assembleMatrixPattern(localAssembler, patternBuilder);
  }

  template <class LA, class MB>
  void assembleMatrixEntries (LA& localAssembler, MB& matrixBackend)
  {
    assembler_.assembleMatrixEntries(localAssembler, matrixBackend);
  }

  template <class LA, class VB>
  void assembleVectorEntries (LA& localAssembler, VB& rhsBackend)
  {
    assembler_.assembleVectorEntries(localAssembler, rhsBackend);
  }

private:
  ColoredRange coloredRange_;
  Executor executor_;
  Assembler assembler_;
};

} // end namespace Dune::BGN

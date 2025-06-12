#pragma once

#include <vector>
#include <string>
#include <tuple>
#include <cmath>
#include <iostream>

#include <fmt/format.h>

#if USE_FLOAT128

template <>
struct fmt::formatter<Dune::Float128>
     : fmt::formatter<double>
{
  template <class FormatContext>
  auto format(Dune::Float128 v, FormatContext& ctx)
  {
    return fmt::formatter<double>::format(double(v), ctx);
  }
};

#endif

namespace Dune {
namespace Impl {

  template <class T, class F1, class F2>
  void printErrorsImpl(std::vector<T> const& widths,
                  std::vector<std::string> const& names,
                  std::vector<std::vector<T>> const& errors,
                  F1 print_line, F2 print_break)
  {
    using std::log;
    // calculate experimental order of convergence
    std::vector<std::vector<T>> eocs(errors.size(), std::vector<T>(1,T(0)));
    for (std::size_t j = 0; j < errors.size(); ++j) {
      for (std::size_t i = 1; i < widths.size(); ++i) {
        eocs[j].push_back( log(errors[j][i]/errors[j][i-1]) / log(widths[i]/widths[i-1]) );
      }
    }

    // print the headlines
    std::vector<std::pair<std::string,std::string>> headlines;
    for (std::size_t j = 0; j < names.size(); ++j)
      headlines.push_back(std::make_pair(names[j], "(eoc) "));

    print_line("lev","h", headlines);
    print_break(headlines.size());

    // print the data lines
    for (std::size_t i = 0; i < widths.size(); ++i) {
      std::vector<std::pair<std::string,std::string>> data;
      for (std::size_t j = 0; j < errors.size(); ++j) {
        data.push_back(std::make_pair(
          fmt::format("{:8.4e}",errors[j][i]),
          fmt::format("{:6.3f}",eocs[j][i]))
        );
      }
      print_line(i, fmt::format("{:<8.6f}",widths[i]), data);
    }
  }

} // end namespace Impl

template <class T>
void printErrorsClassic (std::ostream& out,
                         std::vector<T> const& widths,
                         std::vector<std::string> const& names,
                         std::vector<std::vector<T>> const& errors)
{
  auto print_line = [&out](auto arg0, auto arg1, auto const& args) {
    out << fmt::format("{:<3} | {:<8} ",arg0,arg1);
    for (auto const& arg : args) {
      out << fmt::format("| {:<13} | {:<5} ", arg.first, arg.second);
    }
    out << std::endl;
  };

  auto print_break = [&out](std::size_t n) {
    out << fmt::format("{:-<3}---{:-<8}-","-","-");
    for (std::size_t i = 0; i < n; ++i) {
      out << fmt::format("-{:-<13}---{:-<7}-", "-", "-");
    }
    out << std::endl;
  };

  out << std::endl;
  Impl::printErrorsImpl(widths, names, errors, print_line, print_break);
}


template <class T>
void printErrorsLatex (std::ostream& out,
                       std::vector<T> const& widths,
                       std::vector<std::string> const& names,
                       std::vector<std::vector<T>> const& errors)
{
  auto print_line = [&out](auto arg0, auto arg1, auto const& args) {
    out << fmt::format("{:<3} & {:<8} ",arg0,arg1);
    for (auto const& arg : args) {
      out << fmt::format("& {:<13} & {:<5} ", arg.first, arg.second);
    }
    out << "\\\\" << std::endl;
  };

  auto print_break = [&out](std::size_t n) {
    out << "\\hline" << std::endl;
  };

  out << std::endl;
  out << "\\begin{tabular}{c|c";
  for (std::size_t i = 0; i < names.size(); ++i)
    out << "||c|c";
  out << "}" << std::endl;
  Impl::printErrorsImpl(widths, names, errors, print_line, print_break);
  out << "\\end{tabular}" << std::endl;
}


template <class T>
void printErrors (std::ostream& out,
                  std::vector<T> const& widths,
                  std::vector<std::string> const& names,
                  std::vector<std::vector<T>> const& errors)
{
  auto print_line = [&out](auto arg0, auto arg1, auto const& args) {
    out << fmt::format("{:<3} \u2502 {:<8} ",arg0,arg1);
    for (auto const& arg : args) {
      out << fmt::format("\u2502\u2502 {:<13} \u2502 {:<5} ", arg.first, arg.second);
    }
    out << std::endl;
  };

  auto print_break = [&out](std::size_t n) {
    out << fmt::format("{:\u2500<3}\u2500\u253c\u2500{:\u2500<8}\u2500","\u2500","\u2500");
    for (std::size_t i = 0; i < n; ++i) {
      out << fmt::format("\u2524\u251c\u2500{:\u2500<13}\u2500\u253c\u2500{:\u2500<5}\u2500", "\u2500", "\u2500");
    }
    out << std::endl;
  };

  out << std::endl;
  Impl::printErrorsImpl(widths, names, errors, print_line, print_break);
}

template <class T>
void printErrorsCSV (std::ostream& out,
                     std::vector<T> const& widths,
                     std::vector<std::string> const& names,
                     std::vector<std::vector<T>> const& errors)
{
  auto print_line = [&out](auto arg0, auto arg1, auto const& args) {
    out << arg0 << ", " << arg1;
    for (auto const& arg : args)
      out << ", " << arg.first << ", " << arg.second;
    out << std::endl;
  };

  auto print_break = [&out](std::size_t n) {};

  Impl::printErrorsImpl(widths, names, errors, print_line, print_break);
}

} // end namespace Dune

#pragma once

#include <vector>
#include <string>
#include <tuple>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <stdexcept>

namespace Dune::BGN {

// from https://stackoverflow.com/a/26221725
template<typename ... Args>
std::string string_format( const std::string& format, Args ... args )
{
  int size_s = std::snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
  if( size_s <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
  auto size = static_cast<size_t>( size_s );
  std::unique_ptr<char[]> buf( new char[ size ] );
  std::snprintf( buf.get(), size, format.c_str(), args ... );
  return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}

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
      headlines.push_back(std::make_pair<std::string,std::string>(names[j], "(eoc) "));

    print_line("lev", "h", headlines);
    print_break(headlines.size());

    // print the data lines
    for (std::size_t i = 0; i < widths.size(); ++i) {
      std::vector<std::pair<std::string,std::string>> data;
      for (std::size_t j = 0; j < errors.size(); ++j) {
        std::string
        data.push_back(std::make_pair(
          string_format("%8.4e",errors[j][i]),
          string_format("%6.3f",eocs[j][i]))
        );
      }
      print_line(std::to_string(i), string_format("%8.6f",widths[i]), data);
    }
  }

} // end namespace Impl

template <class T>
void printErrorsClassic (std::ostream& out,
                         std::vector<T> const& widths,
                         std::vector<std::string> const& names,
                         std::vector<std::vector<T>> const& errors)
{
  auto print_line = [&out](std::string arg0, std::string arg1, auto const& args) {
    out << string_format("%3s | %8s ",arg0.c_str(),arg1.c_str());
    for (auto const& arg : args) {
      out << string_format("| %3s | %5s ", arg.first.c_str(), arg.second.c_str());
    }
    out << std::endl;
  };

  auto print_break = [&out](std::size_t n) {
    out << string_format("%s---%s-",std::string(3,'-').c_str(),std::string(8,'-').c_str());
    for (std::size_t i = 0; i < n; ++i) {
      out << string_format("-%s---%s-", std::string(13,'-').c_str(), std::string(7,'-').c_str());
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

// Copyright (C) 2014 Benjamin Kehlet
//
// This file is part of mshr.
//
// mshr is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// mshr is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with mshr.  If not, see <http://www.gnu.org/licenses/>.
//

#include <mshr/STLFileReader.h>

#include <dolfin/geometry/Point.h>
#include <dolfin/common/constants.h>
#include <dolfin/log/LogStream.h>
#include <dolfin/log/log.h>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <string>
#include <map>

namespace
{
  inline double strToDouble(const std::string& s)
  {
    std::istringstream is(s);
    double val;
    is >> val;

    return val;
  }

// A custom compare function class to be used in std::map
// to avoid duplicated vertices.
class FuzzyPointLess
{
 public:
  FuzzyPointLess(double arg_ = 1e-10) : epsilon(arg_) {}
  bool operator()(const std::array<double, 3>& left,
                  const std::array<double, 3>& right) const
  {
    if (std::abs(left[0]-right[0]) > epsilon)
      return left[0] < right[0];

    if (std::abs(left[1] - right[1]) > epsilon)
      return left[1] < right[1];

    if (std::abs(left[2] - right[2]) > epsilon)
      return left[2] < right[2];

    return false;
  }
 private:
  double epsilon;
};

inline void get_next_line(std::ifstream& file, std::string& line, std::size_t &lineno)
{
  std::getline(file, line);
  lineno++;
}

} // end anonymous namespace
//-----------------------------------------------------------------------------
namespace mshr
{

void STLFileReader::read(const std::string filename,
                         std::vector<std::array<double, 3> >& vertices,
                         std::vector<std::vector<std::size_t> >& facets)
{

  dolfin::log(dolfin:: TRACE, "Reading surface from %s ", filename.c_str());

  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

  std::ifstream file(filename.c_str());
  if (!file.is_open())
  {
    dolfin::dolfin_error("PolyhedronUtils.cpp",
                         "open .stl file to read 3D surface",
                         "Failed to open file");
  }

  std::size_t num_vertices = 0;
  std::map<std::array<double, 3>, std::size_t, FuzzyPointLess> vertex_map;
  std::string line;
  std::size_t lineno = 0;
  const boost::char_separator<char> sep(" ");

  // Read the first line and trim away whitespaces
  get_next_line(file, line, lineno);
  boost::algorithm::trim(line);

  if (line.substr(0, 5) != "solid")
  {
    dolfin::dolfin_error("STLFileReader.cpp",
                         "open .stl file to read 3D surface",
                         "File does not start with \"solid\" (line %u", lineno);
  }

  // TODO: Read name of solid  
  get_next_line(file, line, lineno);
  boost::algorithm::trim(line);

  while (file.good())
  {
    bool has_normal = false;
    dolfin::Point normal;

    // Read the line "facet normal n1 n2 n3"
    {
      tokenizer tokens(line, sep);
      tokenizer::iterator tok_iter = tokens.begin();

      if (*tok_iter != "facet")
        dolfin::dolfin_error("STLFileReader.cpp",
                             "open .stl file to read 3D surface",
                             "Expected keyword \"facet\" (line %u)", lineno);
      ++tok_iter;

      // Check if a normal different from zero is given
      if (tok_iter != tokens.end())
      {
        if  (*tok_iter != "normal")
          dolfin::dolfin_error("STLFileReader.cpp",
                               "open .stl file to read 3D surface",
                               "Expected keyword \"normal\"(line %u)", lineno);
        ++tok_iter;

        //dolfin::cout << "Read line: " << line << dolfin::endl;
        
        for (std::size_t i = 0; i < 3; ++i)
        {
          normal[i] = strToDouble(*tok_iter);
          ++tok_iter;
        }

        if (normal.norm() > DOLFIN_EPS)
          has_normal = true;
        
        if (tok_iter != tokens.end())
          dolfin::dolfin_error("PolyhedronUtils.cpp",
                               "open .stl file to read 3D surface",
                               "Expected end of line (line %u)", lineno);
      }
    }

    // if (has_normal)
    // {
    //   dolfin::cout << "Has normal" << dolfin::endl;
    //   dolfin::cout << normal << dolfin::endl;
    // }

    // Read "outer loop" line
    get_next_line(file, line, lineno);
    boost::algorithm::trim(line);

    if (line != "outer loop")
      dolfin::dolfin_error("PolyhedronUtils.cpp",
                           "open .stl file to read 3D surface",
                           "Expected keyword 'outer loop' (line %u)", lineno);

    std::vector<std::size_t> v_indices;

    get_next_line(file, line, lineno);
    boost::algorithm::trim(line);

    tokenizer tokens(line, sep);
    tokenizer::iterator tok_iter = tokens.begin();

    if (*tok_iter != "vertex")
    {
      dolfin::dolfin_error("PolyhedronUtils.cpp",
                           "open .stl file to read 3D surface",
                           "Expected keyword vertex (line %u)", lineno);
    }

    int counter = 0;

    // Read lines with vertices
    do
    {
      // Advance to next
      ++tok_iter;

      const double x = strToDouble(*tok_iter); ++tok_iter;
      const double y = strToDouble(*tok_iter); ++tok_iter;
      const double z = strToDouble(*tok_iter); ++tok_iter;

      std::array<double, 3> vertex = {x, y, z};

      // TODO: Use std::map::find()
      // (to avoid two queries)
      if (vertex_map.count(vertex) > 0)
        v_indices.push_back(vertex_map[vertex]);
      else
      {
        vertex_map[vertex] = num_vertices;
        v_indices.push_back(num_vertices);
        vertices.push_back(vertex);
        num_vertices++;
      }

      // Get next line
      get_next_line(file, line, lineno);
      boost::algorithm::trim(line);

      tokens = tokenizer(line, sep);
      tok_iter = tokens.begin();

      counter++;
    } while (*tok_iter == "vertex");

    // Read 'endloop' line
    boost::algorithm::trim(line);
    if (line != "endloop")
    {
      dolfin::dolfin_error("STLFileReader.cpp",
                           "open .stl file to read 3D surface",
                           "Expected keyword endloop (line %u)", lineno);
    }

    get_next_line(file, line, lineno);
    boost::algorithm::trim(line);
    if (line != "endfacet")
    {
      dolfin::dolfin_error("STLFileReader.cpp",
                           "open .stl file to read 3D surface",
                           "Expected keyword endfacet (line %u)", lineno);
    }

    // Add facet to output
    facets.push_back(v_indices);

    // Get orientation right if normal is given
    if (has_normal) 
    {
      // Compute normal
      const dolfin::Point v1(3, vertices[v_indices[0]].data());
      const dolfin::Point v2(3, vertices[v_indices[1]].data());
      const dolfin::Point v3(3, vertices[v_indices[2]].data());

      const dolfin::Point a = v2-v1;
      const dolfin::Point b = v3-v1;

      dolfin::Point n = a.cross(b);
      n /= n.norm();

      // dolfin::cout << "Normal: " << n << dolfin::endl;
      // if ( (n - normal).norm() > 1e-5 )
      //   dolfin::cout << "Diff: " << (n - normal).norm() << dolfin::endl;
    }

    // Get next line 
    // either start of next facet or endsolid
    get_next_line(file, line, lineno);
    boost::algorithm::trim(line);

    if (line.substr(0, 5) != "facet")
      break;
  }

  // Read the 'endsolid' line
  tokenizer tokens(line, sep);
  tokenizer::iterator tok_iter = tokens.begin();

  if (*tok_iter != "endsolid")
  {
    dolfin::dolfin_error("PolyhedronUtils.cpp",
                         "open .stl file to read 3D surface",
                         "Expected keyword endsolid at line %u", lineno);
  }
  ++tok_iter;

  // TODO: Check name of solid

  dolfin::log(dolfin::TRACE, "Done reading surface");
}

}

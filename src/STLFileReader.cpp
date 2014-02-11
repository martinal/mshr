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
}
//-----------------------------------------------------------------------------
namespace mshr
{

void STLFileReader::read(const std::string filename,
                         std::vector<std::array<double, 3> > vertices,
                         std::vector<std::array<std::size_t, 3> > facets)
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
  std::map<std::array<double, 3>, std::size_t> vertex_map;
  std::string line;
  const boost::char_separator<char> sep(" ");

  // Read the first line and trim away whitespaces
  std::getline(file, line);
  boost::algorithm::trim(line);

  if (line.substr(0, 5) != "solid")
  {
    dolfin::dolfin_error("STLFileReader.cpp",
                         "open .stl file to read 3D surface",
                         "File does not start with \"solid\"");
  }

  // TODO: Read name of solid  
  std::getline(file, line);
  boost::algorithm::trim(line);

  while (file.good())
  {
    //bool has_normal = false;
    //Point normal;

    // Read the line "facet normal n1 n2 n3"
    {
      tokenizer tokens(line, sep);
      tokenizer::iterator tok_iter = tokens.begin();

      if (*tok_iter != "facet")
        dolfin::dolfin_error("STLFileReader.cpp",
                             "open .stl file to read 3D surface",
                             "Expected keyword \"facet\"");
      ++tok_iter;

      // Check if a normal different from zero is given
      if (tok_iter != tokens.end())
      {
        //dolfin::cout << "Expecting normal" << dolfin::endl;

        if  (*tok_iter != "normal")
          dolfin::dolfin_error("STLFileReader.cpp",
                               "open .stl file to read 3D surface",
                               "Expected keyword \"normal\"");
        ++tok_iter;

        //dolfin::cout << "Read line: " << line << dolfin::endl;
        
        // for (std::size_t i = 0; i < 3; ++i)
        // {
        //   normal[i] = strToDouble(*tok_iter);
        //   ++tok_iter;
        // }


        //dolfin::cout << "Normal: " << normal << dolfin::endl;
        // if (normal.norm() > DOLFIN_EPS)
        //   has_normal = true;
        
        // if (tok_iter != tokens.end())
        //   dolfin::dolfin_error("PolyhedronUtils.cpp",
        //                "open .stl file to read 3D surface",
        //                "Expected end of line");
      }
    }

    // Read "outer loop" line
    std::getline(file, line);
    boost::algorithm::trim(line);

    if (line != "outer loop")
      dolfin::dolfin_error("PolyhedronUtils.cpp",
                           "open .stl file to read 3D surface",
                           "Expected key word 'outer loop'");

    std::array<std::size_t, 3> v_indices;

    // Read lines with vertices
    for (std::size_t i = 0; i < 3; ++i)
    {
      std::getline(file, line);
      boost::algorithm::trim(line);

      //dolfin::cout << "read line: " << line << dolfin::endl;

      tokenizer tokens(line, sep);
      tokenizer::iterator tok_iter = tokens.begin();

      if (*tok_iter != "vertex")
      {
        dolfin::dolfin_error("PolyhedronUtils.cpp",
                             "open .stl file to read 3D surface",
                             "Expected key word vertex");
      }
      ++tok_iter;

      const double x = strToDouble(*tok_iter); ++tok_iter;
      const double y = strToDouble(*tok_iter); ++tok_iter;
      const double z = strToDouble(*tok_iter); ++tok_iter;

      std::array<double, 3> vertex = {x, y, z};

      // TODO: Use std::map::find()
      // (to avoid two queries)
      if (vertex_map.count(vertex) > 0)
        v_indices[i] = vertex_map[vertex];
      else
      {
        vertex_map[vertex] = num_vertices;
        v_indices[i] = num_vertices;
        vertices.push_back(vertex);
        num_vertices++;
      }
    }

    // Read 'endloop' line
    std::getline(file, line);
    boost::algorithm::trim(line);
    if (line != "endloop")
    {
      dolfin::dolfin_error("STLFileReader.cpp",
                           "open .stl file to read 3D surface",
                           "Expected key word endloop");
    }

    std::getline(file, line);
    boost::algorithm::trim(line);
    if (line != "endfacet")
    {
      dolfin::dolfin_error("STLFileReader.cpp",
                           "open .stl file to read 3D surface",
                           "Expected key word endfacet");
    }

    // Add facet to output
    facets.push_back(v_indices);

    // Get next line 
    // either start of next facet or endsolid
    std::getline(file, line);
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
                         "Expected key word endsolid");
  }
  ++tok_iter;

  // TODO: Check name of solid

  dolfin::log(dolfin::TRACE, "Done reading surface");
}

}

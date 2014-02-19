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

#include <mshr.h>
#include <dolfin.h>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include <string>
#include <iostream>


namespace po = boost::program_options;

// This program reads in a surface from file and generates mesh
int main(int argc, char** argv)
{
  // Command line options
  po::options_description desc("Generate mesh from surface file");
  desc.add_options()
    ("infile",     po::value<std::string>(), "Filename of infile")
    ("outfile",    po::value<std::string>()->default_value("outmesh.xml"), "Filename of outfile")
    ("resolution", po::value<double>()->default_value(5.0), "Resolution of result mesh")
    ("help",                                                "write help message");

  //Command line positional arguments
  po::positional_options_description positional_options; 
  positional_options.add("infile", 1); 

  // parse command line
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv)
            .options(desc)
            .positional(positional_options).run(), vm);

  // Print help message if requested 
  // (before notify to avoid error messages if only --help is given)
  if (vm.count("help")) { std::cout << desc << std::endl;  exit(0); }

  po::notify(vm);

  std::cout << "Infile: " << vm["infile"].as<std::string>() << std::endl;
  mshr::Surface3D surf(vm["infile"].as<std::string>());

  dolfin::Mesh m;

  mshr::CSGMeshGenerator::generate(m, surf, 30);
  
  dolfin::File f(vm["outfile"].as<std::string>());
  f << m;
    

  std::cout << "Done" << std::endl;

  return 0;
}

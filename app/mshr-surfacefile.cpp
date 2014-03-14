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


namespace
{
void print_mesh_statistics(const dolfin::Mesh& m)
{
  std::cout << "Dolfin mesh of topological dimension " << m.topology().dim() << std::endl;
  std::cout << "  " << m.num_vertices() << " vertices" << std::endl;
  std::cout << "  " << m.num_cells() << " cells" << std::endl;

  const std::pair<double, double> radii_ratio = dolfin::MeshQuality::radius_ratio_min_max(m);
  std::cout << "Minimum cell radii ratio: " << radii_ratio.first  << std::endl;
  std::cout << "Maximum cell radii ratio: " << radii_ratio.second << std::endl;
}
//-----------------------------------------------------------------------------
// Define options and parse command line
void handle_commandline(int argc, char** argv, po::variables_map &vm)
{
  // Command line options
  po::options_description visible("Generate mesh from surface file");
  visible.add_options()
    ("outfile,o",    po::value<std::string>(), "Filename of generated Dolfin mesh")
    ("resolution,r", po::value<double>()->default_value(15.0), "Resolution of result mesh")
    ("stats,s", "Write some statistics of the mesh to stdout")
    ("polyout", po::value<std::string>(), "Write the polyhedron to .poly which Tetgen can read (and do not create a mesh)")
    ("help,h",   "write help message");

  // Options not shown to the user
  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", po::value<std::string>(), "Input file")
    ;

  po::options_description cmdline("Generate mesh from surface file");
  cmdline.add(visible).add(hidden);

  //Command line positional arguments
  po::positional_options_description p;
  p.add("input-file", 1);

  // parse command line
  po::store(po::command_line_parser(argc, argv)
            .options(cmdline)
            .positional(p).run(), vm);

  // Print help message if requested 
  // (before notify to avoid error messages if only --help is given)
  if (vm.count("help"))
  {
    std::cout << visible << std::endl;
    exit(EXIT_SUCCESS);
  }

  po::notify(vm);

}
} //end anonymous namespace
//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  // Ensure Dolfin initializes MPI correctly.
  #ifdef HAS_MPI
  dolfin::SubSystemsManager::init_mpi();
  #endif

  po::variables_map vm;
  handle_commandline(argc, argv, vm);

  // Read the infile
  mshr::Surface3D surf(vm["input-file"].as<std::string>());

  if (vm.count("polyout"))
  {
    // Write the polyhedron to tetgen's file format and exit
    mshr::CSGCGALDomain3D domain(surf);
    mshr::TetgenFileWriter::write(domain,
                            vm["polyout"].as<std::string>());
    exit(EXIT_SUCCESS);
  }


  // Generate the mesh
  dolfin::Mesh m;
  mshr::CSGMeshGenerator::generate(m, surf, vm["resolution"].as<double>());

  // Output mesh if requested
  if (vm.count("outfile"))
  {
    dolfin::File f(vm["outfile"].as<std::string>());
    f << m;
  }
    
  if (vm.count("stats"))
    print_mesh_statistics(m);


  return EXIT_SUCCESS;
}

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
#include <boost/filesystem.hpp>

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

  const std::pair<double, double> volume_min_max = mshr::DolfinMeshUtils::cell_volume_min_max(m);
  std::cout << "Min cell volume: " << volume_min_max.first << std::endl;
  std::cout << "Max cell volume: " << volume_min_max.second << std::endl;
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
    ("polystats", "Write statistics of polyhedron (and do not create a mesh")
    ("backend,b", po::value<std::string>()->default_value("cgal"), "Use 3D mesh generation backend [tetgen|cgal]")
    ("degenerate_tolerance", po::value<double>()->default_value(1e-12), "Tolerance for considering a facet as degenerate. Set to 0 to not remove degenerate facets")
    ("check-mesh", "Check consistency of output mesh (most for debugging/testing")
    ("verbose,v", "Output more information about what is going on")
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

  if (vm.count("verbose"))
    dolfin::set_log_level(dolfin::TRACE);

  // Read the infile
  if (!boost::filesystem::exists(vm["input-file"].as<std::string>()))
  {
    std::cerr << "File " << vm["input-file"].as<std::string>() << "does not exist" << std::endl;
    exit(1);
  }


  mshr::Surface3D surf(vm["input-file"].as<std::string>());
  surf.degenerate_tolerance = vm["degenerate_tolerance"].as<double>();

  // Operations that disable mesh generation
  if (vm.count("polyout") || vm.count("polystats"))
  {
    mshr::CSGCGALDomain3D domain(surf);

    if (vm.count("polyout"))
    {
      std::string extension = boost::filesystem::extension(vm["polyout"].as<std::string>());

      if (extension == ".poly")
      {
        // Write the polyhedron to tetgen's file format
        mshr::TetgenFileWriter::write(domain,
                                      vm["polyout"].as<std::string>());
      }
      else if (extension == ".off")
      {
        domain.save_off(vm["polyout"].as<std::string>());
      }
      else
      {
        std::cerr << "Unknown file type: " << extension << std::endl;
        exit(1);
      }
    }

    if (vm.count("polystats"))
      std::cout << domain.str(true) << std::endl;

    exit(EXIT_SUCCESS);
  }

  // Generate the mesh
  dolfin::Mesh m;

  mshr::generate(m,
                 surf,
                 vm["resolution"].as<double>(),
                 vm["backend"].as<std::string>());

  // Output mesh if requested
  if (vm.count("outfile"))
  {
    dolfin::File f(vm["outfile"].as<std::string>());
    f << m;
  }
    
  if (vm.count("stats"))
    print_mesh_statistics(m);

  if (vm.count("check-mesh"))
  {
    if (!mshr::DolfinMeshUtils::check_mesh(m))
    {
      std::cout << "  Error: Mesh check failed" << std::endl;
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}

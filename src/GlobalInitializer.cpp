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

#include <mshr/GlobalInitializer.h>

#include <CGAL/Random.h>

#include <iostream>

GlobalInitializer::GlobalInitializer()
{
  std::cout << "Initializing globals" << std::endl;
  #ifdef INIT_RANDOM_GENERATOR
    std::cout << "Using fixed random generator generator (seed=" << INIT_RANDOM_GENERATOR << ")" << std::endl;
    CGAL::default_random = CGAL::Random( INIT_RANDOM_GENERATOR );
  #else
    std::cout << "Not fixing seed" << std::endl;
  #endif

}
//-----------------------------------------------------------------------------
GlobalInitializer::~GlobalInitializer()
{

}
//-----------------------------------------------------------------------------
GlobalInitializer _the_singleton;
GlobalInitializer& GlobalInitializer::instance()
{
  return _the_singleton;
}

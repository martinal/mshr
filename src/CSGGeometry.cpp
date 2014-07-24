// Copyright (C) 2012 Anders Logg
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
// Modified by Benjamin Kehlet, 2013


#include <dolfin/common/NoDeleter.h>
#include <mshr/CSGGeometry.h>

namespace mshr
{

//-----------------------------------------------------------------------------
CSGGeometry::CSGGeometry()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
CSGGeometry::~CSGGeometry()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void CSGGeometry::set_subdomain(std::size_t i, std::shared_ptr<CSGGeometry> s)
{
  if (dim() != 2)
    dolfin::dolfin_error("CSGGeometry.cpp",
                 "setting subdomain",
                 "Subdomains are currently supported only in 2D");

  if (s->dim() != dim())
    dolfin::dolfin_error("CSGGeometry.cpp",
                 "setting subdomain",
                 "Subdomain and domain must be of same dimension. Domain was dimension %d and subdomain was %d", dim(), s->dim());

  if (i == 0)
  {
    dolfin::dolfin_error("CSGGeometry.cpp",
                         "Setting reserved CSG subdomain (0)",
                         " Subdomain 0 is reserved and cannot be set by user");
  }

  // Check if i already used
  std::list<std::pair<std::size_t, std::shared_ptr<const CSGGeometry> > >::iterator it = subdomains.begin();
  while (it != subdomains.end())
  {
    if (it->first == i)
    {
      dolfin::warning("Double declaration of CSG subdomain with index %u.", i);

       // Remove existing declaration
       it = subdomains.erase(it);
    }
    else
      ++it;
  }

  subdomains.push_back(std::make_pair(i, s));
}
//-----------------------------------------------------------------------------
void CSGGeometry::set_subdomain(std::size_t i, CSGGeometry& s)
{
  set_subdomain(i, reference_to_no_delete_pointer(s));
}
//-----------------------------------------------------------------------------
bool CSGGeometry::has_subdomains() const
{
  return subdomains.size() > 0;
}

}

====
mshr
====

mshr is a mesh generation component of `FEniCS <http://fenicsproject.org/>`_. It generates tetrahedral `Dolfin <https://bitbucket.org/fenics-project/dolfin>`_ meshes in 2D and 3D (utilizing CGAL and Tetgen as backends) from CSG geometries or surface files.

Installation
============

To build mshr, run::

  cmake <path to mshr source tree>
  make
  make install

mshr's build script will also build CGAL and Tetgen from source and include them in the binary.


Depencies
=========

mshr needs `Dolfin <https://bitbucket.org/fenics-project/dolfin>`_ with python support (pyDolfin). `CGAL <http://www.cgal.org/>`_ and `Tetgen <http://www.tetgen.org>`_ are shipped with mshr and built from source automatically. CGAL needs `Gnu GMP <https://gmplib.org/>`_ and `Gnu MPFR <http://www.mpfr.org/>`_.


License
=======

mshr is licensed under GPL version 3 or (at your option) any later version.


Contact
=======

mshr is hosted at https://bitbucket.org/benjamik/mshr/

For comments and requests, send an email to the FEniCS mailing list::

  fenics@fenicsproject.org

For bug reports and feature requests, visit mshr's issue tracker at BitBucket::

 https://bitbucket.org/benjamik/mshr/issues

Contribute
==========

Contributions are welcome! 

Please read about contributing to FEniCS here: http://fenicsproject.org/contributing/ 

If you plan to implement a new feature, please discuss it at the FEniCS mailing list beforehand. Smaller patches and bugfixes are easiest submitted as `pull request on Bitbucket <https://confluence.atlassian.com/display/BITBUCKET/Work+with+pull+requests>`_.
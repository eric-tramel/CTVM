# CTVM
Implementation of total variation minimization (TVM) in C++. Based on the TVAL3 algorithm of Chengbo Li, "An efficient augmented Lagrangian method with applications to total variation minimization," 2009.

# Dependencies
The Boost libraries, specifically, the Boost uBLAS implementation. 
[See the Boost Website](http://www.boost.org) for details. 

# Building
To build the libraries and test executables, first ensure that you have the following environment variables defined to point to your local Boost installation,

```bash
    BOOST_LIBDIR=/path/to/boost/libs
    BOOST_INCDIR=/path/to/boost/inc
```

For example, in the case of a Homebrew Boost installation on OS X, one
might utilize the common default paths

```bash
    BOOST_LIBDIR=/usr/local/lib
    BOOST_INCDIR=/usr/local/include
```
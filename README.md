# CTVM
Implementation of total variation minimization (TVM) in C++. Based on the TVAL3 algorithm of Chengbo Li, "An efficient augmented Lagrangian method with applications to total variation minimization," 2009.

# Build Dependencies
- *Boost:* The Boost libraries, specifically, the Boost uBLAS implementation. 
[See the Boost Website](http://www.boost.org) for details. 

- *ImageMagick:* ImageMagick [is available here](http://www.imagemagick.org). Specifically, one needs to be able to access the Magick++ header files and associated libraries. 

# Building

To build the libraries and test executables, first ensure that you have the following environment variables defined to point to your local Boost installation,

```bash
    BOOST_LIBDIR=/path/to/boost/libs
    BOOST_INCDIR=/path/to/boost/inc
```
## Unix/Linux/MacOSX
For example, in the case of a Homebrew Boost installation on OS X, one
might utilize the common default paths

```bash
    BOOST_LIBDIR=/usr/local/lib
    BOOST_INCDIR=/usr/local/include
```

Additionally, the Makefile makes use of the `Magick++-config` script distributed with ImageMagick, so one needs to ensure that this script is on the path.

## Windows
It should also be possible to utilize IDE platforms such as Microsoft Visual C++ to compile the project, however, this is untested. The project should be configured such that the ImageMagick and Boost headers and libraries are visible to the project. Otherwise, there is no OS-specific code in this project.

Here, the paths are build in the case of a Cygwin development.

*ImageMagick*
One could simply install particular Cygwin repositories for using directly ImageMagick from the command line:
Graphics > ImageMagick: Image manipulation software suite (utilities)
		 + ImageMagick-doc: ~ ~ ~ ~ (documentation)
		 + libImageMagick1: ~ ~ ~ ~ (runtime)


*Boost*
One easy way to link the Boost headers and libraries in the bash is to create these files in the Cygwin directory (which usually is in the Hard Drive C:). This directory corresponds to the HOME path in the Cygwin command line. Then the bash variable is the same as in the Unix environment

```bash
	BOOST_LIBDIR=/usr/local/lib
	BOOST_INCDIR=/usr/local/include
```
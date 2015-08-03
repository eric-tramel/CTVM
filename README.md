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

### Visual Studio Configuration
Here, we detail the installation and configuration of dependencies on a native Windows environment. 

#### 1. Install and Configure MS VS Community 2015
First, one needs to install Microsoft Visual Studio Community 2015. The IDE can be downloaded
[from here](https://www.visualstudio.com/downloads/download-visual-studio-vs). 
> _Make sure to select the Visual C++ tools as well as the 3rd party command-line tools when performing the installation!_

#### 2. Install ImageMagick
Download the latest Windows version of ImageMagick [from here](http://www.imagemagick.org/script/binary-releases.php#windows). 
Since we will be targeting 64-bit deployment, be sure to download the Win64 version at 16 bits-per-pixel. This should
be the default proposed latest ImageMagick release.

> _When installing ImageMagick, be sure to choose the option to install C++ header files and development tools!_

For ease, choose an easy install path. For these instructions we choose `C:\ImageMagick` as the install path for 
the ImageMagick tools.

#### 3. Install and Build Boost
Download the latest relase of Boost [from here](http://sourceforge.net/projects/boost/files/boost/1.58.0/). 
We have currently tested that the build works for Boost 1.58.0. 

1. After the download completes, extract the compressed Boost download to any directory.
2. Next, open `Developer Command Prompt for VS2015` and navigate the directory Boost was extracted to.
3. Execute the command `bootstrap.bat` to build the Boost configuration utilities.
4. Decide on an install directory for Boost. For ease we will choose `C:\boost` for these instructions.
5. Execute the following command to build and install Boost in 64 bit mode `.\b2 install --prefix=C:\boost toolset=msvc address-model=64`
6. The install process may take some time, so one must be patient ! 

#### 4. Set Environment Variables
In order to allow for minimal changes to the Visual Studio solution and project files, we refer to the 
library and include search paths for both Boost and Imagemagick using environment variables. This
means that prior to running the solution for the first time, one needs to set the following
environment variables:
```bash
BOOST_INCDIR = C:\path\to\boost\include
BOOST_LIBDIR = C:\path\to\boost\lib
IMAGEMAGICK_INCDIR = C:\path\to\imagemagick\include
IMAGEMAGICK_LIBDIR = C:\path\to\imagemagick\lib
```
This can be accomplished by going to
`Control Panel > System and Security > System` and then selecting
"Advanced System Settings" and pressing the "Environment Variable"
button at the bottom.
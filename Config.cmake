############################
# Specify the c++ compiler. icpx and g++ have been tested.

## For linux system:
## Uncomment the following line if you use an g++ compiler
set(ENV{CC} gcc)
set(ENV{CXX} g++)


## For MacOS, do not use gcc or g++, because they are alias of Apple clang and clang++ instead of the real
## GNU-gcc. In addition, apple clang may cause some problems in compiling, e.g. missing openmp. I'd suggest using
## homebrew clang.
## To use the GNU gcc and g++, first intall them with homebrew, and use the command gcc-[version] or g++-[version]
# set(ENV{CC} clang)
# set(ENV{CXX} clang++)

# Uncomment the following two lines if you use an intel compiler
#set(ENV{CXX} icpx)
#set(USE_MKL TRUE)


############################netcdf support###############################
# Specify whether netcdf library will be used. The value is either TRUE or FALSE
# Uncomment the following line if you want to write results to netcdf files. Comment it if you don't use netcdf.
set(USE_NETCDF TRUE)


##############################location of the final binary files###########
#set (CMAKE_INSTALL_PREFIX "prefix_path")

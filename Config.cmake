############################
# Specify the c++ compiler. icpx and g++ have been tested.
# Uncomment the following line if you use an g++ compiler
set(ENV{CXX} g++)

# Uncomment the following two lines if you use an intel compiler
#set(ENV{CXX} icpx)
#set(USE_MKL TRUE)


############################netcdf support###############################
# Specify whether netcdf library will be used. The value is either TRUE or FALSE
# Uncomment the following line if you want to write results to netcdf files. Comment it if you don't use netcdf.
set(USE_NETCDF TRUE)


##############################location of the final binary files###########
#set (CMAKE_INSTALL_PREFIX "prefix_path")

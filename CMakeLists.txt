# Set the minimum version of CMake that can be used
# To find the cmake version run
# $ cmake --versionh
cmake_minimum_required(VERSION 3.10)

# set the project name
project(GraAda3D)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

find_package(OpenMP REQUIRED)

message("${CMAKE_CXX_COMPILER_ID}")
message("Compiler command: $ENV{CXX}")
message("Openmp flag ${OpenMP_CXX_FLAGS}")
if (DEFINED ENV{USE_NETCDF})
    message("USE_NETCDF=$ENV{USE_NETCDF}; The netcdf library will be included.")
endif()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
#message("${CMAKE_CXX_COMPILER}")
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message("Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# Forward modelling code
add_library(Forward_obj OBJECT
    src/Forward_modelling/Cell.cpp
    src/Forward_modelling/Fwd.cpp
    src/Forward_modelling/Face.cpp
    src/Forward_modelling/GravFormula.cpp
    src/Forward_modelling/Integral.cpp
    src/Forward_modelling/Mesh.cpp
    src/Forward_modelling/Observation.cpp
    src/Forward_modelling/Point.cpp
    src/Forward_modelling/RectPrism.cpp)
target_include_directories(Forward_obj
    PUBLIC
        src/Forward_modelling
        3rd_party_lib/boost_1_72_0
        3rd_party_lib/eigen
        3rd_party_lib/Linterp
        3rd_party_lib/Timer
        #if netcdf library is used
        $<$<BOOL:$ENV{USE_NETCDF}>:
            3rd_party_lib/netcdf/netcdf-cxx4-4.3.1
            3rd_party_lib/netcdf/netcdf-c-4.7.4
        >
        )

target_link_libraries(Forward_obj
        PUBLIC
        OpenMP::OpenMP_CXX
        $<$<BOOL:$ENV{USE_NETCDF}>:netcdf_c++4>
        gsl
        gslcblas        
        )
    
target_compile_options(Forward_obj
    PUBLIC
    -Wfatal-errors 
    ${OpenMP_CXX_FLAGS}
    )
target_compile_definitions(Forward_obj
    PUBLIC
    $<$<CXX_COMPILER_ID:Intel>:USE_MKL>
    $<$<BOOL:$ENV{USE_NETCDF}>:USE_NETCDF>
    )
# inversion code
add_library(Inversion_obj OBJECT
    src/Inversion/InversionBase.cpp
    src/Inversion/AdaptiveInversion.cpp
    src/Inversion/LinearInversion.cpp
    src/Inversion/GaussNewtonInversion.cpp)
target_include_directories(Inversion_obj
    PUBLIC
        src/Inversion
        )
target_link_libraries(Inversion_obj
        PUBLIC 
        Forward_obj
        )
############################################################
# Create an executable
############################################################
add_executable(GraAda3D
    src/GraAda3D.cpp
)
target_include_directories(GraAda3D
    PRIVATE
        src/Forward_modelling
        src/Inversion
        )

target_link_options(GraAda3D
    PRIVATE
    $<$<BOOL:$ENV{USE_MKL}>:-L$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl>
    )
target_link_libraries(GraAda3D
    PRIVATE
    Forward_obj
    Inversion_obj
    )
add_executable(Test_ordering
    src/Test_ordering.cpp)
target_link_libraries(Test_ordering
    PRIVATE
    Forward_obj
    Inversion_obj
    )

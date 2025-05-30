# GraAda3D

## 1 Introduction

It is a 3D gravity inversion program implemented with C++. The inversion domain is discretized as rectangular prismatic meshes in the Cartesian coordinate system. The inversion mesh can be adaptively refined to boost computational performance and avoid over-parameterization. In addition, users can choose to use L0-norm (minimum support) and L1-norm regularization to obtain a focused image, or use L2-norm regularization to get a smooth result. Furthermore, a-priori information can be incorporated in inversion through cross-gradient coupling or direct parameter relationship. 

Any combination of gravity field components or gravity gradient components (gz, gx, gy, Tzz, Txz, Tyz, Txx, Txy, Tyy) can be used as input data.  Exact analytical solutions of gravity field and gravity gradient tensor are used to ensure accuracy of the forward modeling. 

**Hightlights**

- Adaptively refined mesh for inversion parameters
- Incorporation of  a-priori constraints through cross-gradient coupling or direct parameter relation
- Lp-norm regularization (p=0,1,2,...)

## 2 Installation

### 2.1 Dependencies

(1) [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page): A C++ template library for linear algebra. A copy of Eigen is contained in the `./3rd_party_lib/` folder and you don't need to install it manually.

(2) [Linterp](https://rncarpio.github.io/linterp/): A multidimensional interpolation class which is dependent on the Boost library. A copy of Linterp is contained in `./3rd_party_lib` folder so that you don't need to install it manually.

(3) [netcdf](https://www.unidata.ucar.edu/software/netcdf/) [*optional*]: A library to write and read data in .nc format. NetCDF is a popular format to store and share data. Data of this format can be read in GMT (generic mapping tool) or python. The installation of netcdf is optional. If you don't need netcdf, you can edit `GraAda3D/Config.cmake` and change  `set(USE_NETCDF TRUE)` to `set(USE_NETCDF FALSE)`.

The fastest way to install netcdf is using the built-in package manager:

```shell
#For ubuntu users
sudo apt-get install libnetcdf-dev
sudo apt-get install libnetcdf-c++4-dev

#For Fedora users
sudo yum install netcdf-devel
sudo yum install netcdf-cxx-devel
```

> For CentOS user, it is recommended to build netcdf-c and netcdf-cxx from the source code. See https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html for netcdf-c and https://github.com/Unidata/netcdf-cxx4 for the netcdf c++ library.

(4) GNU Scientific Library (gsl): **[gsl](https://www.gnu.org/software/gsl/)** is used to calculate wavelet transforms. `cd` the directory containing the source code of gsl, follow the instructions in file INSTALL (`./configure`, `make`, `make check`, `make install`, `make installcheck`).

### 2.2 Buiding

#### Build with CMake (recommended)

Before compiling the code, open file `Config.cmake`, follow the instruction in the comments to edit it.

```
# go to the GraAda3D folder
cd GraAda3D

# remove the "build" directory if it exists
rm -rf ./build

# Create an empty directory named build
mkdir build

cd build

# Generate a makefile
cmake ..

# Compile
make

# copy the executable files to GraAda3D/bin
make install
```

Then you can add GraAda3D/bin to the environmental variable PATH. First, go to the GraAda3D/bin directory and type `pwd`, which will show you the absolute path of the our bin directory:

```shell
cd ../bin

pwd -P
```
For example, the `pwd` command on my computer shows `/home/zhongyiyuan/GraAda3D/bin`. Then add the following line to file `~/.bashrc`
```
export PATH=/home/zhongyiyuan/GraAda3D/bin${PATH:+:${PATH}}
```

Now start a new terminal, type `GraAda3D` to check whether it has been installed successfully.

#### Build with Make

Build the program using Make,

```shell
cd GraAda3D
make
```

By default, the program is built without netcdf support. In this way, you don't have to install the netcdf library in advance, but there are not *.nc files generated after inversion. The results are written in  *.vtk file, which can be visulaized by [Paraview](https://www.paraview.org/).


The default compiler is g++. To build with an [Intel C++ compiler](https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit/download.html) (`icpc`, `icpx`, `dpcpp`), use the CXX option,

```shell
make CXX=icpx
```

> The latest Intel compiler is Intel oneAPI for which the command is `icpx` or `dpcpp`. The command of a classic intel compiler is `icpc`. It seems that the program compiled with the intel compiler has a slightly better performance than g++.


## 3 Examples

### 3.1 Forward modelling

We have implemented some handy tool for gravity forward modelling, including:

- `makeModel`: design a synthetic model.
- `generatePoints`: generate a list of observation points
- `compute`: model the gravity and/or gravity gradient of a synthetic model at given observation points. A file specifying the model and a file specifying the observation points are required. The model file can be generated by `makeModel`. The observation points can be generated by `generatePoints`.

An example is given in `Examples/Forward_modelling/model1`.

### 3.2 Inversion

Now, you have 2 executable programs in the GraAda3D folder: **GraAda3D** and **Synthetic_data1**.

The **GraAda3D** program is the  inversion program, and **Synthetic_data1** is used to generate a synthetic data set for validation of the inversion program. 

Change from the current directory to `GraAda3D/Examples` and type the following command in the terminal:

```
../Synthetic_data1
```

to generate synthetic data with noise. Then we get a text file named *dobs_g_z* which contains vertical gravity data $g_z$, and 3 files named *dobs_T_zz*, *dobs_T_xz*, *dobs_T_yz* which contains gravity gradient components Tzz, Txz, Tyz. 

### 3.2.1 A example of inversion using gz data

1. cd `GraAda3D/Examples/Inv_gz` In this folder, we have prepared the following configuration files:  **config**, **config_data**, **config_inversion** and **config_model**. 
   
   The file **config** is the only parameter required by the inversion program, which specifies another 3 configuration files: **config_data**, **config_inversion** and **config_model**. 
   
   - The file **config_data** contains information about data. 
   - The file **config_inversion** specifies controlling parameters for the inversion. 
   - The file **config_model** specifies  the inversion region (size and initial mesh discretization). 
   
   These configuration files allow comments that start with '#'.  See the comments in the configuration files for detailed explanation. 

2. Then, let's run the inversion

```
../../GraAda3D config
```

After the inversion is finished,  the result is written into **gz_result.vtk** file. The vtk file can be visualized in [Paraview](https://www.paraview.org/). 

If `GraAda3D` is compiled with NetCDF (`make USE_NETCDF=1`), the resulting model is also written into a ***.nc** file. A python script `crossSections2.py`  is used to read the result in *.nc file and plot cross sections.

3. Similarly, `cd GraAda3D/Examples/Inv_Tzz` and run the inversion of Tzz data with`../../GraAda3D config`.

### 3.2.2 Inversion using  multiple components of gravity gradient tensor

The file **config_data** tells the program which data to be used. 

1. Enter `Inv_Txz_Tyz_Tzz` folder,  open **config_data**, and compare it with the ones in `Inv_Tzz` folder and `Inv_gz` folder.
2. Run the inversion

```
../../GraAda3D config
```

3. Show the inversion result in Paraview, or using python:

```
python crossSections2.py
python Txz_Tyz_Tzz.py
```

### 3.2.3 Examples of inversion with a priori information

A priori information (e.g. existing velocity models) can be included in the inversion to improve the reliability of the inversion result.  The a priori model should be gridded data given in a XYZ table.

#### (1) Reference and initial model

Using direct empirical relations between density and the known parameter, we can obtain a reference density model. In this algorithm, if a reference model is specified, it will also serve as the initial model. 

We have prepared a file **ref_model** in the GraAda3D/Examples folder, which contains gridded data of a reference model.

1. cd `GraAda3D/Examples/Inv_with_ref`, open the file **config_inversion**, look at line 70, which specify the reference model and see comments on lines 63-69 for explanation.

2. Run

```
../../GraAda3D config
```

#### (2) Cross gradient constraint

1. cd `GraAda3D/Examples/Inv_with_cross_gradient_constraint`, open **config_inversion**, look at line 65. 

2. Run

```
../../GraRect config
```

### 3.2.4 Inversion using L1-norm regularization

1. Enter `L1_inv_gz` folder,  open **config_inversion**, and compare it with the one in `Inv_gz` folder. See Line 14 in **config_inversion**, it contains 2 parameters about Lp-norm inversion. The first parameter is the "p" of Lp-norm, the second parameter is a small value which is used to avoid singularity. In the file `L1_inv_gz/config_inversion`, line 14, p=1, which means L1-norm regularization is used here.
2. Run the inversion

```
../../GraAda3D config
```

3. Show the inversion result in Paraview, or using python:

```
python crossSections2.py
```

### 3.2.5 Inversion using minimum support regularization

1. Enter `L0_inv_gz(Minimum_support)`, open **config_inversion**, see line 14. Here, p=0.

2. Run the inversion,

```
../../GraAda3D config
```



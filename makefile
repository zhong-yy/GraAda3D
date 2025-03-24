CXX:= g++
CXX:=$(strip $(CXX))
ifeq ($(CXX),g++)
CXXFLAGS:=-std=c++11 -O3 -fopenmp -Wfatal-errors #-g
endif
ifeq ($(CXX),icpc)
CXXFLAGS:=-std=c++11 -O3 -qopenmp -Wfatal-errors -mkl=parallel
endif
ifeq ($(CXX),icpx)
CXXFLAGS:=-std=c++11 -O3 -qopenmp -Wfatal-errors -I"${MKLROOT}/include" -qmkl=parallel -D USE_MKL    
endif
ifeq ($(CXX),dpcpp)
CXXFLAGS:=-std=c++11 -O3 -qopenmp -Wfatal-errors -I"${MKLROOT}/include" -mkl=parallel -D USE_MKL    
endif


LINKFLAGS:=
USE_NETCDF:=0

ifeq ($(OS),Windows_NT)
LINKFLAGS += -static
endif
ifeq ($(CXX),icpc)
LINKFLAGS+=-liomp5 -lpthread -lm -ldl
endif
ifeq ($(CXX),icpx)
LINKFLAGS+=  -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
endif
ifeq ($(CXX),dpcpp)
LINKFLAGS+=  -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
endif

LINKFLAGS+= -lgsl -lgslcblas -lm

#libraries
#DIR_EIGEN               :=./3rd_party_lib/Eigen3.3.7
DIR_EIGEN               :=./3rd_party_lib/eigen
DIR_LINTERP             :=./3rd_party_lib/Linterp
DIR_BOOST               :=./3rd_party_lib/boost_1_72_0
INCLUDE_EIGEN           :=-I$(DIR_EIGEN)
INCLUDE_LINTERP         :=-I$(DIR_LINTERP)
INCLUDE_BOOST           :=-I$(DIR_BOOST)
INCLUDE_TIMER           :=-I./3rd_party_lib/Timer
INCLUDE_LIBRARIES         :=-I$(DIR_EIGEN) -I$(DIR_LINTERP) -I$(DIR_BOOST) $(INCLUDE_TIMER)

ifeq ($(USE_NETCDF),1)
CXXFLAGS+=-D USE_NETCDF
LINKFLAGS   +=-lnetcdf_c++4
DIR_NETCDF_CXX          :=./3rd_party_lib/netcdf/netcdf-cxx4-4.3.1
DIR_NETCDF_C            :=./3rd_party_lib/netcdf/netcdf-c-4.7.4
INCLUDE_NETCDF          :=-I$(DIR_NETCDF_CXX) -I$(DIR_NETCDF_C)
INCLUDE_LIBRARIES         +=$(INCLUDE_NETCDF)
endif


#My codes
DIR_INVERSION           :=./src/Inversion
DIR_FWD                 :=./src/Forward_modelling

INCLUDE_SRC             :=-I$(DIR_INVERSION) -I$(DIR_FWD)
SRCS             := $(wildcard $(DIR_INVERSION)/*.cpp)
SRCS             += $(wildcard $(DIR_FWD)/*.cpp)
OBJS             := $(patsubst %.cpp, %.o, $(SRCS))
SRC_PROGRSM_OPTIONS     :=$(wildcard ./3rd_party_lib/boost_1_72_0/libs/program_options/src/*.cpp)
OBJS_PROGRAM_OPTIONS    :=$(patsubst %.cpp, %.o, $(SRC_PROGRSM_OPTIONS))
OBJS             += $(OBJS_PROGRAM_OPTIONS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_LIBRARIES) $(INCLUDE_SRC) $ -c $< -o $@

EXE:=Synthetic_data1 GraAda3D Padding_test makeModel test_wavelet_forward1 test_wavelet_forward2
ALL:$(EXE)
.PHONY: all
Padding_test:$(OBJS) ./src/Padding_test.o
	$(CXX) $(CXXFLAGS) -o Padding_test ./src/Padding_test.o $(OBJS) $(LINKFLAGS)
Synthetic_data1:$(OBJS) ./src/Synthetic_data1.o
	$(CXX) $(CXXFLAGS) -o Synthetic_data1 ./src/Synthetic_data1.o $(OBJS) $(LINKFLAGS)
GraAda3D:$(OBJS) ./src/GraAda3D.o
	$(CXX) $(CXXFLAGS) -o GraAda3D ./src/GraAda3D.o $(OBJS) $(LINKFLAGS)
makeModel:$(OBJS) ./src/makeModel.o
	$(CXX) $(CXXFLAGS) -o makeModel ./src/makeModel.o $(OBJS) $(LINKFLAGS)
test_wavelet_forward1:$(OBJS) ./src/test_wavelet_forward1.o
	$(CXX) $(CXXFLAGS) -o test_wavelet_forward1 ./src/test_wavelet_forward1.o $(OBJS) $(LINKFLAGS)
test_wavelet_forward2:$(OBJS) ./src/test_wavelet_forward2.o
	$(CXX) $(CXXFLAGS) -o test_wavelet_forward2 ./src/test_wavelet_forward2.o $(OBJS) $(LINKFLAGS)

.PHONY: clean
clean:
	@rm -rf *.o *~  $(OBJS) $(EXE) ./src/GraAda3D.o ./src/Synthetic_data1.o ./src/Padding_test.o ./src/makeModel.o ./src/test_wavelet_forward1.o ./src/test_wavelet_forward2.o

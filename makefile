CXX:= g++ 
CXX:=$(strip $(CXX))
ifeq ($(CXX),g++)
CXXFLAGS:=-std=c++11 -O3 -fopenmp -Wfatal-errors #-g
endif
ifeq ($(CXX),icpc)
CXXFLAGS:=-std=c++11 -O3 -qopenmp -Wfatal-errors -mkl=parallel
endif
ifeq ($(CXX),icpx)
CXXFLAGS:=-std=c++11 -O3 -qopenmp -Wfatal-errors -I"${MKLROOT}/include" -mkl=parallel -D USE_MKL    
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

#libraries
#DIR_EIGEN               :=./contrib/Eigen3.3.7
DIR_EIGEN               :=./contrib/eigen
DIR_LINTERP             :=./contrib/Linterp
DIR_BOOST               :=./contrib/boost_1_72_0
INCLUDE_EIGEN           :=-I$(DIR_EIGEN)
INCLUDE_LINTERP         :=-I$(DIR_LINTERP)
INCLUDE_BOOST           :=-I$(DIR_BOOST)
INCLUDE_TIMER           :=-I./contrib/Timer
INCLUDE_CONTRIB         :=-I$(DIR_EIGEN) -I$(DIR_LINTERP) -I$(DIR_BOOST) $(INCLUDE_TIMER)

ifeq ($(USE_NETCDF),1)
CXXFLAGS+=-D USE_NETCDF
LINKFLAGS   +=-lnetcdf_c++4
DIR_NETCDF_CXX          :=./contrib/netcdf/netcdf-cxx4-4.3.1
DIR_NETCDF_C            :=./contrib/netcdf/netcdf-c-4.7.4
INCLUDE_NETCDF          :=-I$(DIR_NETCDF_CXX) -I$(DIR_NETCDF_C)
INCLUDE_CONTRIB         +=$(INCLUDE_NETCDF)
endif


#My codes
DIR_INVERSION           :=./src/Inversion
DIR_FWD                 :=./src/Forward_modelling

INCLUDE_SRC             :=-I$(DIR_INVERSION) -I$(DIR_FWD)
SRCS             := $(wildcard $(DIR_INVERSION)/*.cpp)
SRCS             += $(wildcard $(DIR_FWD)/*.cpp)
OBJS             := $(patsubst %.cpp, %.o, $(SRCS))

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_CONTRIB) $(INCLUDE_SRC) $ -c $< -o $@

EXE:=Synthetic_data1 GraInvRect
ALL:$(EXE)
.PHONY: all
Synthetic_data1:$(OBJS) ./src/Synthetic_data1.o
	$(CXX) $(CXXFLAGS) -o Synthetic_data1 ./src/Synthetic_data1.o $(OBJS) $(LINKFLAGS)
GraInvRect:$(OBJS) ./src/GraInvRect.o
	$(CXX) $(CXXFLAGS) -o GraInvRect ./src/GraInvRect.o $(OBJS) $(LINKFLAGS)

.PHONY: clean
clean:
	@rm -rf *.o *~  $(OBJS) $(OBJS_TESSEROID) $(EXE) ./src/GraInvRect.o ./src/Synthetic_data1.o

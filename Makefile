PREFIX=/Users/amedhi/Projects/Codes/SimpleVMC
# Libraries
EIGEN_INCLUDE=-I/usr/local/include
EIGEN_USE_MKL=EIGEN_USE_MKL_ALL
MKL_INCLUDE=-I/opt/intel/mkl/include/intel64/lp64
MKL_LDFLAGS=-L/opt/intel/mkl/lib
MKL_LIBS=-lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
INCLUDE = $(EIGEN_INCLUDE)

# Compiler
CXX=g++ -std=c++11 
CPPFLAGS= #-D$(EIGEN_USE_MKL)
OPTFLAGS=-Wall -O3
CXXFLAGS=$(CPPFLAGS) $(OPTFLAGS) $(INCLUDE)
LDFLAGS=$(MKL_LDFLAGS) 
LIBS=$(MKL_LIBS)

BUILD_DIR=$(PREFIX)/build

#-------------------------------------------------------------
# Source files
SRC = main.cpp
SRC+= lattice.cpp
SRC+= random.cpp
SRC+= basis.cpp
SRC+= wavefunction.cpp
SRC+= sysconfig.cpp
SRC+= mcdata/mcdata.cpp
SRC+= mcdata/mc_observable.cpp
SRC+= vmc.cpp
SRCS=$(addprefix src/,$(SRC))
#-------------------------------------------------------------
# Headers
HDR = constants.h
HDR+= lattice.h
HDR+= random.h
HDR+= basis.h
HDR+= matrix.h
HDR+= wavefunction.h
HDR+= sysconfig.h
HDR+= mcdata/mcdata.h
HDR+= mcdata/mc_observable.h
HDR+= vmc.h
HDRS=$(addprefix src/,$(HDR))
#-------------------------------------------------------------
# Target
TAGT=a.out

# All .o files go to BULD_DIR
OBJS=$(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SRCS))
# GCC/Clang will create these .d files containing dependencies.
DEPS=$(patsubst %.o,%.d,$(OBJS)) 

.PHONY: all
all: $(TAGT) #$(INCL_HDRS)

$(TAGT): $(OBJS)
	$(CXX) -o $(TAGT) $(OBJS) $(LDFLAGS) $(LIBS)  

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -o $@ $<

# Include all .d files
-include $(DEPS)

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	@echo "$(CXX) -c $(CXXFLAGS) -o $(@F) $(<F)"
	@$(CXX) -MMD -c $(CXXFLAGS) -o $@ $<

.PHONY: clean
clean:	
	@echo "Removing temporary files in the build directory"
	@rm -f $(OBJS) $(DEPS) 
	@echo "Removing $(TAGT)"
	@rm -f $(TAGT) 


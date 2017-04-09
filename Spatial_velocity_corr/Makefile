
# Location of home
HOME=/usr/users/iff_th2/duman

# Location of the scripts
S_PATH=$(HOME)/Cells_in_LAMMPS/Scripts 
VPATH=../Utility

# Location of the HDF5 package
H5_PATH=$(HOME)/hdf5_parallel
H5_INCLUDE=$(H5_PATH)/include/
H5_LIB=$(H5_PATH)/lib/

# Compiler and delete options
CXX=g++
RM=rm -f

# Compilfer flags
DEFAULT_FLAGS=-O3 
PARFLAGS=-fopenmp
LDFLAGS=-I$(H5_INCLUDE) -L$(H5_LIB)
LDLIBS=-lhdf5 -lm
CXXFLAGS=$(DEFAULT_FLAGS) $(LDFLAGS) $(LDLIBS)

# Sources and objects
SRCS=calc_sp_velocity_corr.cpp read_write.cpp
OBJS=calc_sp_vel_corr.o read_write.o

### note that "$<" gets the first prerequisite
### $@ is the target
# Main target
calc_sp_vel_corr: $(OBJS)
	$(CXX) $(CXXFLAGS) -o calc_sp_vel_corr $(OBJS) 

# Targets
calc_sp_vel_corr.o: calc_sp_velocity_corr.cpp read_write.hpp basic.hpp
		 $(CXX) -c $(CXXFLAGS) $< -o $@

read_write.o: read_write.cpp basic.hpp
	      $(CXX) -c $(CXXFLAGS) $< -o $@

# Cleaning options
.PHONY: clean
clean:
	$(RM) $(OBJS)
	$(RM) calc_sp_vel_corr


CXX = g++
MPICXX = mpic++

#source code in another directory
SRC = ../src
vpath %.cc $(SRC) #vpath usage: vpath <pattern> <directories>
CXXFLAGS = -Wall -I../src
OMPFLAG  = -fopenmp -lpthread

pde_solver_obj = grid.o pde_solver.o driver_pde_solver.o
pde_solver_omp_obj =  grid.o pde_solver.o pde_solver_omp.o driver_pde_solver_omp.o
test_grid_obj  = grid.o test_grid.o 
pde_solver_ompi_src = driver_pde_solver_ompi.cc pde_solver_ompi.cc grid.cc grid_ompi.cc
pde_solver_ompi_src := $(addprefix $(SRC)/, $(pde_solver_ompi_src))
exe_list = heat_serial test_grid.x heat_omp heat_mpi

all: $(exe_list) 

heat_serial: $(pde_solver_obj)
	$(CXX) -o $@ $^

test_grid.x: $(test_grid_obj)
	$(CXX) -o $@ $^

heat_omp: $(pde_solver_omp_obj)
	$(CXX) $(OMPFLAG) -o $@ $^

heat_mpi: $(pde_solver_ompi_src)
	$(MPICXX) -o $@ $^

%.o: %.cc
	$(CXX) $(OMPFLAG) $(CXXFLAGS) -c $< 

clean:
	$(RM) *.o *.x


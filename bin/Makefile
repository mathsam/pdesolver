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
exe_list = pde_solver.x test_grid.x pde_solver_omp.x pde_solver_ompi.x

all: $(exe_list) 

pde_solver.x: $(pde_solver_obj)
	$(CXX) -o $@ $^

test_grid.x: $(test_grid_obj)
	$(CXX) -o $@ $^

pde_solver_omp.x: $(pde_solver_omp_obj)
	$(CXX) $(OMPFLAG) -o $@ $^

pde_solver_ompi.x: $(pde_solver_ompi_src)
	$(MPICXX) -o $@ $^

%.o: %.cc
	$(CXX) $(OMPFLAG) $(CXXFLAGS) -c $< 

clean:
	$(RM) *.o *.x


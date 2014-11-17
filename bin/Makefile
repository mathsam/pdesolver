CXX = g++

#source code in another directory
SRC = ../src
vpath %.cc $(SRC) #vpath usage: vpath <pattern> <directories>
CXXFLAGS = -Wall -I../src

pde_solver_obj = grid.o driver_pde_solver.o pde_solver.o
test_grid_obj  = grid.o test_grid.o 
exe_list = pde_solver.x test_grid.x

all: $(exe_list) 

pde_solver.x: $(pde_solver_obj)
	$(CXX) -o $@ $^

test_grid.x: $(test_grid_obj)
	$(CXX) -o $@ $^

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< 

clean:
	$(RM) *.o *.x


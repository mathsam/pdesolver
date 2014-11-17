#include"pde_solver_omp.h"
#include<iostream>
#include<cstdlib>

static const double kPI = 3.14159265359;

int main(int argc, char* argv[]){
    if(argc != 3){
        std::cout << "USAGE: "<< argv[0]<< " <nx> <num_thread>"<< std::endl;
        std::exit(EXIT_FAILURE);
    }

    int nx = std::atoi(argv[1]);
    int num_thread = std::atoi(argv[2]);
    int ny = nx;
    double dx = kPI / double(nx);
    double dy = kPI / double(ny);
    double dt = dx*dx / 8.0;

    double stopping_time = kPI*kPI / 2.0;

    PdeSolverOmp mysolver(nx,ny,dt,dx,dy,num_thread);
    mysolver.RunSimulation(stopping_time);
}

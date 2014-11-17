#include"pde_solver.h"
#include<iostream>
#include<cstdlib>

static const double kPI = 3.14159265359;

int main(int argc, char* argv[]){
    if(argc != 2){
        std::cout << "USAGE: "<< argv[0]<< " <nx>"<< std::endl;
        std::exit(EXIT_FAILURE);
    }

    int nx = std::atoi(argv[1]);
    int ny = nx;
    double dx = kPI / double(nx);
    double dy = kPI / double(ny);
    double dt = dx*dx / 8.0;

    double stopping_time = kPI*kPI / 2.0;

    PdeSolver mysolver(nx,ny,dt,dx,dy);
    mysolver.RunSimulation(stopping_time);
}

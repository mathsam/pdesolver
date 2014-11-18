#include"pde_solver_ompi.h"
#include<iostream>
#include<cstdlib>
#include"mpi.h"

static const double kPI = 3.14159265359;

int main(int argc, char* argv[]){
    int rank, num_proc;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(argc != 2){
        std::cout << "USAGE: "<< argv[0]<< " <nx>"<< std::endl;
        std::exit(EXIT_FAILURE);
    }

    int nx = std::atoi(argv[1]);
    if(nx%num_proc != 0){
        std::cout << "Domain not divisible by number of processors"  << std::endl;
        std::exit(EXIT_FAILURE);
    }
    int nx_per_proc = nx/num_proc;

    int ny = nx;
    double dx = kPI / double(nx);
    double dy = kPI / double(ny);
    double dt = dx*dx / 8.0;
    double stopping_time = kPI*kPI / 2.0;
    int halo_width_x = 1;
    int root = 0;

    PdeSolverOmpi mysolver(nx_per_proc, ny, rank, num_proc, halo_width_x,
                           dt,dx,dy);
    mysolver.RunSimulation(stopping_time, root);
}

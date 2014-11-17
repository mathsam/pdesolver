#include"pde_solver.h"

static const double kPI = 3.14159265359;

int main(){
    int nx = 200;
    int ny = 200;
    double dx = kPI / double(nx);
    double dy = kPI / double(ny);
    double dt = dx*dx / 8.0;

    double stopping_time = kPI*kPI / 2.0;

    PdeSolver mysolver(nx,ny,dt,dx,dy);
    mysolver.RunSimulation(stopping_time);
}

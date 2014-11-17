#ifndef _PDE_SOLVER_OMP_H
#define _PDE_SOLVER_OMP_H
#include"pde_solver.h"

/**
 * @brief a framework to solve 2d PDE in parallel using OpenMPI
 *
 * Inherited from serial class PdeSolver 
 */

class PdeSolverOmp: public PdeSolver{
public:
    /**
     * @breif generator for solver
     * @param nx number of grid points in x direction
     * @param ny number of grid points in y direction
     * @param dt time step for time-stepping
     * @param dx grid spacing in x direction
     * @param dy grid spacing in y direction
     * @param num_thread number of thread to parallelize
     */
    PdeSolverOmp(int nx, int ny, double dt, double dx, double dy, int num_thread);

    void TimeStepping();

private:
    const int num_thread_;
};

#endif //_PDE_SOLVER_OMP_H

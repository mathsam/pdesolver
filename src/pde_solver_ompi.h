#ifndef _PDE_SOLVER_OMPI_H
#define _PDE_SOLVER_OMPI_H
#include"grid_ompi.h"
/**
 * @brief a framework to solve 2d PDE using Open MPI to parallel 
 *
 * Mainly the same as Grid for serial program, but uses GridOmpi
 * Main tasks are 
 * (1) spacial differentiation 
 * (2) time stepping
 * (3) I/O
 */

class PdeSolverOmpi{
public:
    /**
     * @breif generator for solver
     * @param nx number of grid points in x direction for each processor
     * @note nx is not the size of global grid
     * @param ny number of grid points in y direction
     * @param rank current processor rank
     * @param num_proc total number of processors
     * @param halo_width_x width of halo grid in x direction
     * @param dt time step for time-stepping
     * @param dx grid spacing in x direction
     * @param dy grid spacing in y direction
    */
    PdeSolverOmpi(int nx, int ny, int rank, int num_proc, int halo_width_x,
                  double dt, double dx, double dy);

    virtual void RunSimulation(double stopping_time, int root=0);

    ///spacial derivative at each grid point for time stepping
    double DfDt(int ix, int jy);

    ///time stepping
    virtual void TimeStepping();

    ///write field to disk
    void WriteField(int root=0);
protected:
    const double dt_;
    const double dx_;
    const double dy_;
    const int nx_wo_halo_; ///<number of grid without halo grid
    const int ny_;
    const int rank_;
    const int num_proc_;
    double current_time_;
    GridOmpi field2d_;

    ///domain to be time stepped
    int min_x_;
    int max_x_;
    int min_y_;
    int max_y_;
};

#endif //_PDE_SOLVER_OMPI_H

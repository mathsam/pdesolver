#ifndef _PDE_SOLVER_H
#define _PDE_SOLVER_H
#include"grid.h"
/**
 * @brief a framework to solve 2d PDE 
 *
 * Main tasks are 
 * (1) spacial differentiation 
 * (2) time stepping
 * (3) I/O
 */

class PdeSolver{
public:
    /**
     * @breif generator for solver
     * @param nx number of grid points in x direction
     * @param ny number of grid points in y direction
     * @param dt time step for time-stepping
     * @param dx grid spacing in x direction
     * @param dy grid spacing in y direction
    */
    PdeSolver(int nx, int ny, double dt, double dx, double dy);

    virtual void RunSimulation(double stopping_time);

    /**
     *@brief set up the initial condition
     *
     * here simply set everything to zero
    */
    virtual void InitializeField();

    ///spacial derivative at each grid point for time stepping
    double DfDt(int ix, int jy);

    ///time stepping
    virtual void TimeStepping();

    ///write field to disk
    void WriteField();
protected:
    const int nx_;
    const int ny_;
    const double dt_;
    const double dx_;
    const double dy_;
    double current_time_;
    Grid field2d_;

    ///domain to be time stepped
    int min_x_;
    int max_x_;
    int min_y_;
    int max_y_;
};

#endif //_PDE_SOLVER_H

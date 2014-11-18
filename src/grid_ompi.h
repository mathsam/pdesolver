#ifndef _SMART_GRID_OMPI_H_
#define _SMART_GRID_OMPI_H_
#include"grid.h"
/**
 * @brief smart grid handles openmpi calls to update its halo under
 * the hood, and again handles the boundary condition.
 *
 * Each object has a rank_, which corresponds to the Open MPI rank.
 * Only handles domain decomposition in the x-direction at this time.
 */

class GridOmpi: public Grid{
public:
    /**
     *@brief domain size nx x ny plus halo grid on the left and right
     * with size halo_width_x
     * 
     * The inner domain with size nx x ny is what time stepping is working
     * on; while halo grid is from adjacent domain.
     *@param nx number of grids in x
     *@param ny number of grids in y
     *@param rank current rank of processor
     *@param num_proc total number of processors
     *@param halo_width_x size of halo, default is 1
    */
    GridOmpi(int nx, int ny, int rank, int num_proc,
             int halo_width_x = 1);

    ///only part of the domain of the original problem. so need to override
    ///the base class
    void InitializeBoundary();

    ///because one does not need to solve for halo, so min_x starts from 1 now
    void DomainToSolve(int & min_x, int & max_x,
                       int & min_y, int & max_y);

    void UpdateHalo();

    void GatherField(double* global_field, int root);

private:
    const int rank_;
    const int num_proc_;
    const int halo_width_x_;
    const int left_rank_;
    const int right_rank_;
    const int nx_wo_halo_; ///<size of domain in x direction without halo
}; 

#endif //_SMART_GRID_OMPI_H_

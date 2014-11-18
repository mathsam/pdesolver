#include"grid_ompi.h"
#include<cmath>
#include"mpi.h"

static const double kPI = 3.14159265359;

GridOmpi::GridOmpi(int nx, int ny, int rank, int num_proc,
                   int halo_width_x):
    Grid(nx+2*halo_width_x, ny), 
    rank_(rank), num_proc_(num_proc),
    halo_width_x_(halo_width_x),
    left_rank_(  (rank+num_proc-1)%num_proc),
    right_rank_( (rank+1)%num_proc),
    nx_wo_halo_(nx){
    InitializeBoundary();    
};

void GridOmpi::DomainToSolve(int & min_x, int & max_x,
                   int & min_y, int & max_y){
    min_x = 1; ///< 0 and nx_-1 are halos
    max_x = nx_ - 2;
    min_y = 1; ///< y = 0 and y = ny_-1 are fixed
    max_y = ny_ - 2;
}

void GridOmpi::InitalizeBoundary(){
    for(int ix = 0; ix < nx_; ix++){
        double x = double(rank_*nx_wo_halo_-1 + ix) 
                 / double(nx_wo_halo_ * num_proc_-1) * kPI;
        double sin2_x = std::sin(x) * std::sin(x);
        field2d_A_[ix * ny_]         = 1.0 - sin2_x; ///lower boundary
        field2d_A_[ny_-1 + ix * ny_] = sin2_x;       ///upper boundary
        field2d_B_[ix * ny_]         = 1.0 - sin2_x; ///lower boundary
        field2d_B_[ny_-1 + ix * ny_] = sin2_x;       ///upper boundary
    }
}

/**
 * @brief Update halo grid for the current working field (field that used 
 * by get_point
 *
 * Use UpdateHalo before each time stepping
 */
void GridOmpi::UpdateHalo(){
    int rc;
    MPI_Status Stat;
    ///update all the right halo
    if(rank_%2 == 0){    
       rc = MPI_Send(&field2d_A_[ny_], ny_, MPI_DOUBLE, 
                     left_rank_,///< send to its left domain
                     0, ///<tag for sending field, which is right halo for 
                        ///<its left domain 
                     MPI_COMM_WORLD);
       rc = MPI_Recv(&field2d_A_[(nx_-1)*ny_], ny_, MPI_DOUBLE,
                     right_rank_,/// recieve from its right domain
                     0, MPI_COMM_WORLD, &Stat);
    }
    else{
       rc = MPI_Recv(&field2d_A_[(nx_-1)*ny_], ny_, MPI_DOUBLE,
                     right_rank_,/// recieve from its right domain
                     0, MPI_COMM_WORLD, &Stat); 
       rc = MPI_Send(&field2d_A_[ny_], ny_, MPI_DOUBLE,
                     left_rank_, 0, MPI_COMM_WORLD);
    }

    ///update all the left halo
    if(rank_%2 == 0){    
       rc = MPI_Send(&field2d_A_[(nx_-2)*ny_], ny_, MPI_DOUBLE, 
                     right_rank_, 
                     1, ///<tag for sending field, which is left halo for 
                        ///<its right domain 
                     MPI_COMM_WORLD);
       rc = MPI_Recv(&field2d_A_[0], ny_, MPI_DOUBLE,
                     left_rank_,/// recieve from its left domain
                     1, MPI_COMM_WORLD, &Stat); 
    }
    else{
       rc = MPI_Recv(&field2d_A_[0], ny_, MPI_DOUBLE,
                     left_rank_,/// recieve from its left domain
                     1, MPI_COMM_WORLD, &Stat); 
       rc = MPI_Send(&field2d_A_[(nx_-2)*ny_], ny_, MPI_DOUBLE, 
                     right_rank_, 
                     1, ///<tag for sending field, which is left halo for 
                        ///<its right domain 
                     MPI_COMM_WORLD);
    }

}

/**
 * @brief gather all the domain and put it into one field from all processors
 */
void GridOmpi::GatherField(double* global_field, ///< where to store the gathered field
                           int root///<which processor do the gathering
                          ){
    int rc;
    rc = MPI_Gather(&field2d_A_[ny_], (nx_-2)*ny_, MPI_DOUBLE,
                    global_field, (nx_-2)*ny_, MPI_DOUBLE, root, MPI_COMM_WORLD);
}

#include"pde_solver_omp.h"
#include"omp.h"
#include<stdexcept>
#include<iostream>

PdeSolverOmp::PdeSolverOmp(int nx, int ny, double dt, double dx, double dy, int num_thread):
   PdeSolver(nx, ny, dt, dx, dy), num_thread_(num_thread){
    if(num_thread < 1) 
        throw std::invalid_argument("thread number must be at least 1");
}

void PdeSolverOmp::TimeStepping(){
    if(nx_%num_thread_ != 0)
        throw std::invalid_argument("Domain is not divisible by the number of threads");

    int chunk_columns;
    chunk_columns = nx_/num_thread_;

    omp_set_num_threads(num_thread_);
    ///fork
    #pragma omp parallel default(shared)
    {

        #pragma omp for
        for(int i = 0; i <= max_x_; ++i){
            for(int j = min_y_; j <= max_y_; ++j){
                field2d_.set_point(i,j) = field2d_.get_point(i,j) + DfDt(i,j) * dt_;
            }
        }
    }///merge 
    field2d_.UpdateField();
}


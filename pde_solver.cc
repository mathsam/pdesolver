#include"pde_solver.h"
#include<iostream>
#include<fstream>

static const double kKappa = 1.0;

PdeSolver::PdeSolver(int nx, int ny, double dt, double dx, double dy):
    nx_(nx), ny_(ny),
    dt_(dt), dx_(dx), dy_(dy),
    current_time_(0.0),
    field2d_(nx,ny){
    field2d_.DomainToSolve(min_x_, max_x_,
                           min_y_, max_y_);
};

void PdeSolver::RunSimulation(double stopping_time){
    while(current_time_ < stopping_time){
        TimeStepping();
        current_time_ += dt_;
    }
    WriteField();
}


double PdeSolver::DfDt(int ix, int jy){
    static const double inverse_dxdy = kKappa/dx_/dy_;
    double spacial_derivative;
    spacial_derivative = inverse_dxdy*(field2d_.get_point(ix-1, jy)
                                     + field2d_.get_point(ix+1, jy)
                                     + field2d_.get_point(ix, jy-1)
                                     + field2d_.get_point(ix, jy+1)
                               - 4.0 * field2d_.get_point(ix,   jy));
    return spacial_derivative;
}

void PdeSolver::TimeStepping(){
    for(int i = min_x_; i <= max_x_; ++i){
        for(int j = min_y_; j <= max_y_; ++j){
            field2d_.set_point(i,j) += DfDt(i,j) * dt_; ///simple explicit Euler scheme
        }
    }
}

void PdeSolver::WriteField(){
    std::ofstream result_file;
    result_file.open("result.csv");
    for(int j = 0; j < ny_; ++j){
        for (int i = 0; i < nx_; ++i){
            result_file << field2d_.get_point(i,j) << '\t';
        }
        result_file << '\n';
    }
    result_file.close();
}

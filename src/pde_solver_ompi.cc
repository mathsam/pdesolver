#include"pde_solver_ompi.h"
#include<iostream>
#include<fstream>
#include<iomanip>
#include<iostream>
#include<sstream>
#include<string>
#include<string.h>

static const double kKappa = 1.0;

PdeSolverOmpi::PdeSolverOmpi(int nx, int ny, 
           int rank, int num_proc, int halo_width_x,
           double dt, double dx, double dy):
    dt_(dt), dx_(dx), dy_(dy), 
    nx_wo_halo_(nx), ny_(ny), rank_(rank), num_proc_(num_proc),
    current_time_(0.0),
    field2d_(nx,ny,rank,num_proc,halo_width_x){
    field2d_.DomainToSolve(min_x_, max_x_,
                           min_y_, max_y_);
    std::cout << "Rank: " << rank_ << " Domain: [" << min_x_ 
              << ", " << max_x_ << "]" << std::endl;
};

void PdeSolverOmpi::RunSimulation(double stopping_time, int root){
    InitializeField();
    while(current_time_ < stopping_time){
        field2d_.UpdateHalo();///update halo grid before time stepping
        TimeStepping();
        current_time_ += dt_;
    }
    WriteField(root);
}

void PdeSolverOmpi::InitializeField(){
    for(int i = min_x_; i <= max_x_; ++i){
        for(int j = min_y_; j <= max_y_; ++j){
            field2d_.set_point(i,j) = 0.0; ///simply set the working field to 0
        }
    }
    field2d_.UpdateField();
}

double PdeSolverOmpi::DfDt(int ix, int jy){
    static const double inverse_dxdy = kKappa/dx_/dy_;
    double spacial_derivative;
    spacial_derivative = inverse_dxdy*(field2d_.get_point(ix-1, jy)
                                     + field2d_.get_point(ix+1, jy)
                                     + field2d_.get_point(ix, jy-1)
                                     + field2d_.get_point(ix, jy+1)
                               - 4.0 * field2d_.get_point(ix,   jy));
    return spacial_derivative;
}

void PdeSolverOmpi::TimeStepping(){
    for(int i = min_x_; i <= max_x_; ++i){
        for(int j = min_y_; j <= max_y_; ++j){
            field2d_.set_point(i,j) = field2d_.get_point(i,j) + DfDt(i,j) * dt_; ///simple explicit Euler scheme
        }
    }
    field2d_.UpdateField();
}

void PdeSolverOmpi::WriteField(int root){
    double * global_field = NULL;
    field2d_.UpdateHalo(); //synchronize before output
    if(root == rank_) global_field = new double[nx_wo_halo_ * ny_ * num_proc_];
    field2d_.GatherField(global_field, root);

    if(root == rank_){
        std::ofstream result_file;
        {///formating the output filename
            std::stringstream ss_filename;
            ss_filename << "temp_reso" << nx_wo_halo_ * num_proc_ << "x"
                        << ny_ << "_numproc_" << num_proc_ << ".csv";
            std:: string filename = ss_filename.str();
            char char_filename[filename.size()+1];
            strcpy(char_filename, filename.c_str());
            result_file.open(char_filename);
        }
        result_file << std::setprecision(5) << std::fixed;

        double ave_field = 0;
        for(int j = 0; j < ny_; ++j){
            for (int i = 0; i < nx_wo_halo_ * num_proc_; ++i){
                result_file << global_field[j + i*ny_] << '\t';
                ave_field += global_field[j + i*ny_];
            }
            result_file << '\n';
        }
        result_file.close();
        ave_field = ave_field/double(nx_wo_halo_ * num_proc_ * ny_);
        std::cout << "Average temparature is " << ave_field << std::endl;
    }
}

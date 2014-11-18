#include"grid.h"
#include<stdexcept>
#include<cmath>
#include<new>

static const double kPI = 3.14159265359;


Grid::Grid(int nx, int ny):
    nx_(nx), ny_(ny), field2d_A_(NULL), field2d_B_(NULL){
    if(nx < 2 || ny < 2)
        throw std::invalid_argument("Domain must be at least 2x2");

    field2d_A_ = new double[nx*ny];
    field2d_B_ = new double[nx*ny];
    InitializeBoundary();
}

Grid::Grid(int n): 
    nx_(n), ny_(n), field2d_A_(NULL), field2d_B_(NULL){
    if(n < 2)
        throw std::invalid_argument("Domain must be at least 2x2");

    field2d_A_ = new double[n*n];
    field2d_B_ = new double[n*n];
    InitializeBoundary();
}

Grid::~Grid(){
    delete [] field2d_A_;
    delete [] field2d_B_;
}

void Grid::DomainToSolve(int & min_x, int & max_x,
                   int & min_y, int & max_y){
    min_x = 0;
    max_x = nx_ - 1; 
    min_y = 1; /// y = 0 and y = ny_-1 are fixed
    max_y = ny_ - 2;
}

void Grid::InitializeBoundary(){
    for(int ix = 0; ix < nx_; ix++){
        double x = double(ix) / double(nx_-1) * kPI;
        double sin2_x = std::sin(x) * std::sin(x);
        field2d_A_[ix * ny_]         = 1.0 - sin2_x; ///lower boundary
        field2d_A_[ny_-1 + ix * ny_] = sin2_x;       ///upper boundary
    }
}

double Grid::get_point(int ix, int jy) const{
    Grid* local_this = const_cast<Grid*>(this);
    double point_ij = local_this->GetFieldPoint(field2d_A_, ix, jy);
    return point_ij;
}

double & Grid::set_point(int ix, int jy){
    return GetFieldPoint(field2d_B_, ix, jy);
}

void Grid::UpdateField(){
    double * tmp;
    tmp = field2d_A_;
    field2d_A_ = field2d_B_;
    field2d_B_ = tmp;
}

/**
 * @brief returns the value at grid point (ix, jy)
 * @note spacial rule for negative index or index out of range
 *       in the x direciton due to the periodic boundary condition
*/
double & Grid::GetFieldPoint(double* field2d, int ix, int jy){
    if(ix == 0      || ix == nx_)     return get_left_boundary(field2d, jy);
    if(ix == nx_ -1 || ix == -1)      return get_right_boundary(field2d,jy);
    if(jy == 0      )     return get_lower_boundary(field2d, ix);
    if(jy == ny_ -1 )     return get_upper_boundary(field2d, ix);

    if(ix < 0 || ix >= nx_ || jy < 0 || jy >= ny_)
        throw std::invalid_argument("Grid index outof range");

    return field2d[jy + ix * ny_];
}   

inline double & Grid::get_left_boundary(double* field2d, int jy){
    return field2d[jy];
}

inline double & Grid::get_right_boundary(double* field2d, int jy){
    return field2d[jy + (nx_-1)*ny_];
}

inline double & Grid::get_lower_boundary(double* field2d, int ix){
    return field2d[ix*ny_];
}

inline double & Grid::get_upper_boundary(double* field2d, int ix){
    return field2d[ny_-1 + ix * ny_];
}

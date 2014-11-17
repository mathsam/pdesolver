#ifndef _SMART_GRID_H
#define _SMART_GRID_H

/**
 * @brief grid which stores the field to work on and automatically
 *        handles the boundary condition
 *
 * specifically handles periodic boundary condition in x direction,
 * that is T(0,y) = T(pi,y)
 * and T(x,0) = cos^2(x), T(x,pi) = sin^2(x)
 * indexing follows the C notation that first element has index 0
 */

class Grid{
public:
    ///construct a nx x ny 2 dimensnional field
    Grid(int nx, int ny);

    ///construct a n x n 2 dimensional field
    Grid(int n);

    ~Grid();

    ///set the constant boundary condtion at y = 0 and y = pi
    void InitializeBoundary();

    ///returns the index range for the domain to solve, which is
    /// [min_x, max_x] x [min_y, max_y]
    void DomainToSolve(int & min_x, int & max_x,
                       int & min_y, int & max_y);

    ///returns the value at point (ix, jy)
    double get_point(int ix, int jy) const;

    ///returns the reference for point (ix, jy) so that one can set
    ///the point as set_point(ix,jy) = some_value;
    double & set_point(int ix, int jy);

    inline unsigned int get_dimen_x(){
        return nx_;
    }

    inline unsigned int get_dimen_y(){
        return ny_;
    }

private:
    /**
     * to store 2d field; implemented as an 1d array
     * indexing is faster in y direction, so
     * field2d_(ix, jy) = field2d_[jy + ix * ny_]
    */
    double * field2d_;    
    const int nx_;    ///< number of grids in x direction
    const int ny_;    ///< number of grids in y direction

    double & get_left_boundary (int jy);
    double & get_right_boundary(int jy);
    double & get_lower_boundary(int ix);
    double & get_upper_boundary(int ix);

    ///hide copy constructor and copy assigment operator
    Grid(Grid& );
    Grid& operator=(const Grid& );
};

#endif //_SMART_GRID_H

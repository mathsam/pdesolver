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
 *
 * stores two fields, field_A_ and field_B_
 * One stores current field, the other stores the field one step
 * forward. UpdateField() switches the two field, intending for a
 * time step forward
 */

class Grid{
public:
    ///construct a nx x ny 2 dimensnional field
    Grid(int nx, int ny);

    ///construct a n x n 2 dimensional field
    Grid(int n);

    ~Grid();

    ///set the constant boundary condtion at y = 0 and y = pi
    virtual void InitializeBoundary();

    ///returns the index range for the domain to solve, which is
    /// [min_x, max_x] x [min_y, max_y]
    virtual void DomainToSolve(int & min_x, int & max_x,
                       int & min_y, int & max_y);

    ///returns the value at point (ix, jy)
    double get_point(int ix, int jy) const;

    ///returns the reference for point (ix, jy) so that one can set
    ///the point as set_point(ix,jy) = some_value;
    ///set the field other than used by get_point
    double & set_point(int ix, int jy);

    ///after a time step, update
    void UpdateField();

    inline unsigned int get_dimen_x(){
        return nx_;
    }

    inline unsigned int get_dimen_y(){
        return ny_;
    }

protected:
    const int nx_;    ///< number of grids in x direction
    const int ny_;    ///< number of grids in y direction

    /**
     * to store 2d field; implemented as an 1d array
     * indexing is faster in y direction, so
     * field2d_(ix, jy) = field2d_[jy + ix * ny_]
     * Use two grids to store current field and a time step forward
    */
    double * field2d_A_;    
    double * field2d_B_; 

    double & GetFieldPoint(double* field2d, int ix, int jy);
    double & get_left_boundary (double* field2d, int jy);
    double & get_right_boundary(double* field2d, int jy);
    double & get_lower_boundary(double* field2d, int ix);
    double & get_upper_boundary(double* field2d, int ix);

private:
    ///hide copy constructor and copy assigment operator
    Grid(Grid& );
    Grid& operator=(const Grid& );
};

#endif //_SMART_GRID_H

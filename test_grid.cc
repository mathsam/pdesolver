#include<iostream>
#include<iomanip>
#include"grid.h"

using namespace std;

void print_field(Grid & mygrid){
    cout << setprecision(3);
    for(int j=0; j< mygrid.get_dimen_y(); j++){
        for(int i=0; i< mygrid.get_dimen_x(); i++){
            cout << mygrid.get_point(i,j);
            cout << '\t';
        }
        cout << '\n';
    }
}


int main(){
    int nx = 10;
    int ny = 5;
    Grid mygrid(nx,ny);

    

    cout << "Test boundary condition" <<endl;
    print_field(mygrid);

    cout << "\nTest setting value" << endl;
    int min_x, max_x, min_y, max_y;
    mygrid.DomainToSolve(min_x, max_x, min_y, max_y);
    cout << "\nDomain to solve is xrange x yrange: ["
         <<min_x<<", "<<max_x<<"] x ["
         <<min_y<<", "<<max_y<<"]\n";
    for(int i=min_x; i<=max_x; ++i){
        for(int j=min_y; j<=max_y; ++j){
            mygrid.set_point(i,j) = double(i+j);
        }
    }
    cout << "Setting working field[i,j] = i + j\n";
    print_field(mygrid);
}


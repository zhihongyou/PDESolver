#include <iostream> 
#include <vector>
#include <fstream>
#include <random>
#include <math.h>
#include "fieldclass.h"

using namespace std;


// -----------------------------------------------------------------------
void Field::initFieldConst(Mesh mesh_t, double*& f_t, double f_value) {
    f_t = new double[mesh_t.getGridNumberAllBulk()];
    for (int i=0; i<mesh_t.grid_number.x; i++) {
        for (int j=0; j<mesh_t.grid_number.y; j++) {
            int idx=i*mesh_t.grid_number.x+j;
            f_t[idx]=f_value;
        };
    };
};

// -----------------------------------------------------------------------
void Field::initFieldGaus(Mesh mesh_t, double*& f_t, double r_center, double r_decay, double gaus_amplitude) {
    f_t = new double[mesh_t.getGridNumberAllBulk()];
    for (int i=0; i<mesh_t.grid_number.x; i++) {
        for (int j=0; j<mesh_t.grid_number.y; j++) {
            int idx=i*mesh_t.grid_number.x+j;
            f_t[idx]=gaus_amplitude;
        };
    };
};

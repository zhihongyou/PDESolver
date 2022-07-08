#ifndef FIELDFUNCTION_H
#define FIELDFUNCTION_H
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include <cmath>

using namespace std; 


class Field;
// ---------------------------------------------------------------------
struct FFuncArgs {
    int di;
    int dj;
    int Nx;
    int Ny;
    int Nbx;
    int Nby;
    double dx;
    double dy;
    double f_func_arg1;
    double f_func_arg2;
    double f_func_arg3;

    FFuncArgs(){};

    FFuncArgs (double f_func_arg1_1) {
        f_func_arg1=f_func_arg1_1;
    };

    FFuncArgs (double f_func_arg1_1, double f_func_arg2_1) {
        f_func_arg1=f_func_arg1_1;
        f_func_arg2=f_func_arg2_1;
    };

    FFuncArgs (double f_func_arg1_1, double f_func_arg2_1, double f_func_arg3_1) {
        f_func_arg1=f_func_arg1_1;
        f_func_arg2=f_func_arg2_1;
        f_func_arg3=f_func_arg3_1;
    };
    
    FFuncArgs (int Nx1, int Ny1, int Nbx1, int Nby1, double dx1, double dy1) {
        Nx=Nx1;
        Ny=Ny1;
        Nbx=Nbx1;
        Nby=Nby1;
        di=Nx1+2*Nbx1;
        dj=1;
        dx=dx1;
        dy=dy1;
    };
};

struct FFuncDef {
    // ******************************************************************
    // Struct used to construct a single function call to a field
    // f_perator is a string used to determine function call
    // field_ptr is the field on which the f_operator act.
    // ******************************************************************
    
    string f_operator;
    Field* field_ptr;
    FFuncArgs f_func_args;

    FFuncDef (Field* field_ptr1) {
        f_operator="1";
        field_ptr=field_ptr1;
    };

    FFuncDef (string f_operator1, Field* field_ptr1) {
        f_operator=f_operator1;
        field_ptr=field_ptr1;
    };
    
    FFuncDef (string f_operator1, Field* field_ptr1, FFuncArgs f_func_args1) {
        f_operator=f_operator1;
        field_ptr=field_ptr1;
        f_func_args=f_func_args1;
    };    
};


// -----------------------------------------------------------------------
struct rhsTerm {
    // ******************************************************************
    // Operation struct for collecting a single RHS term of a field.
    // f_funcs provides information on all functions of fields, which are
    //   to be multiplied to form the RHS term.
    // prefactor is the prefactor of this term.
    // scheme determines whether this term is calculated explicitly or
    //   implicitly
    // ******************************************************************
    
    vector<FFuncDef> f_funcs;
    double prefactor;
    string scheme;

    rhsTerm(vector<FFuncDef> f_funcs1) {
        prefactor=1;
        f_funcs=f_funcs1;
        scheme="explicit";
    };

    rhsTerm(double prefactor1, vector<FFuncDef> f_funcs1) {
        prefactor=prefactor1;
        f_funcs=f_funcs1;
        scheme="explicit";
    };
    
    rhsTerm(double prefactor1, vector<FFuncDef> f_funcs1, string scheme1) {
        prefactor=prefactor1;
        f_funcs=f_funcs1;
        scheme=scheme1;
    };
};


// ----------------------------------------------------------------------
struct rhsPtrs {
    // ******************************************************************
    // Struct that determines the RHS of a field
    // num_terms is the number of terms on the RHS
    // prefactors are collections of prefactors of each term.
    // num_func_1term is the number of fields in each term
    // f_func_ptrs collect all pointers on the RHS.
    // ******************************************************************
    int* num_terms;
    double* prefactors;
    int* num_funcs_1term;
    double** f_func_ptrs;
    int* schemes;
};


using FFuncType =  double (*) (double*,int,FFuncArgs);

// Map to determine FFuncType, i.e. which field function, to call based on
//   the operator (first string), and the scheme (second scheme).
map<pair<string,string>, FFuncType> f_func_map_all;
map<pair<string,string>, FFuncType> f_func_map_all_dev;

struct FFuncItem {
    string f_operator;
    int f_func_idx;
    FFuncType f_func;
    FFuncArgs f_func_args;
};


// ----------------------------------------------------------------------
FFuncType getFFuncDevPtr(FFuncType* f_func_ptr_dev) {
    FFuncType f_func_ptr_host;
    cudaMemcpyFromSymbol(&f_func_ptr_host, *f_func_ptr_dev, sizeof(FFuncType));
    return f_func_ptr_host;
};


#endif

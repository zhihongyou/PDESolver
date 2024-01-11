#ifndef FIELDCLASS_EXPORT_CU
#define FIELDCLASS_EXPORT_CU


// -------------------------------------------------------------
void Field::export_conf(string str_t, string device, int include_boun=0) {
    if (device=="cpu") {
        export_conf_any(f[0],name(),str_t, "cpu", include_boun);
    } else if (device=="gpu") {
        export_conf_any(f[0],name(),str_t, "gpu", include_boun);
    };
}


// -------------------------------------------------------------
void Field::export_f_func(string f_operator, string str_t, string device, int include_boun=0) {

    int f_func_idx;
    for (int i=0; i<num_f_funcs; i++) {
        if (f_funcs_rhs[i].f_operator==f_operator) {
            f_func_idx=f_funcs_rhs[i].f_func_idx;
        };
    };
    
    if (device=="cpu") {        
        export_conf_any(f_funcs_host[f_func_idx], name()+"_"+f_operator, str_t, "cpu", include_boun);
    } else if (device=="gpu") {
        export_conf_any(f_funcs_host[f_func_idx], name()+"_"+f_operator, str_t, "gpu", include_boun);
    };
}


// -------------------------------------------------------------
void Field::export_conf_any(double* f_t, string f_name, string str_t, string location_t="cpu" , int include_boun=0) {
    ofstream conf_file;
    int PrecData=8;
    string conf_file_name=dire_expo()+f_name+"_"+ str_t + ".dat";
    conf_file.open(conf_file_name.c_str() );
    
    int idx;
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
    int* idx0=new int [4];

    
    if (location_t=="gpu") {
        allocField<double>(f_temp_host, "cpu");
        cudaMemcpy(f_temp_host, f_t, gridNumberAll()*sizeof(double),cudaMemcpyDeviceToHost);
    };        
    
    if (include_boun==0) {
        idx0[0]=0;
        idx0[1]=0;
        idx0[2]=0;
        idx0[3]=0;
    } else {
        idx0[0]=-Nbx;
        idx0[1]=Nbx;
        idx0[2]=-Nby;
        idx0[3]=Nby;
    };

    for (int j=idx0[2]; j<Ny+idx0[3]; j++) {
        for (int i=idx0[0]; i<Nx+idx0[1]; i++) {        
            idx=(Nx+2*Nbx)*(j+Nby)+i+Nbx;
            if (location_t=="cpu") {
                conf_file<<setiosflags(ios::scientific) <<setprecision(PrecData) <<f_t[idx]<<endl;
            } else {
                conf_file<<setiosflags(ios::scientific) <<setprecision(PrecData) <<f_temp_host[idx]<<endl;
            };
        }
    }
    conf_file.close();
}

// -------------------------------------------------------------

#endif

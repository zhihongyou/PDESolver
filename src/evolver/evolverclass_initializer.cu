#ifndef EVOLVERCLASS_INITIALIZER_CU
#define EVOLVERCLASS_INITIALIZER_CU


// -------------------------------------------------------------
void Evolver::initEvolver() {
    // Initialize field function map used to identify functions to call
    //   based on the operator name.
    setFFuncMap();
    
    // Initialize fields
    initFields();
    
    // Initialize 
    initRHSs();
};


// -------------------------------------------------------------
// Allocate memory for fields 
// Location of memory depends on device of the evolver.
void Evolver::initFields () {

    // Initiate fields and LHS&RHS
    for (auto f_ptr_i : (*system_ptr).field_ptrs ) {
        // Set numerical scheme for each field.
        (*f_ptr_i).FDMScheme=FDMScheme;

        if ((*f_ptr_i).initCond()=="import") {
            int te=floor(time_start/time_export);    
            string str_t=to_string(te);
            // cout << "initial time=" << str_t<<endl;
            (*f_ptr_i).initFieldImport(str_t, 0);
        };
        
        int num_grid=(*f_ptr_i).gridNumberAll();
        if (device=="cpu") {            
            for (int i_f_copy=0; i_f_copy<num_field_copy; i_f_copy++) {
                (*f_ptr_i).f[i_f_copy]=new double[num_grid];
                if ((*f_ptr_i).priority()==0) {
                    (*f_ptr_i).rhs[i_f_copy]=new double[num_grid];
                    (*f_ptr_i).lhs[i_f_copy]=new double[num_grid];
                };
            };
            // Copy data from f_host to device
            for (int idx=0; idx<num_grid; idx++) {
                (*f_ptr_i).f[0][idx]=(*f_ptr_i).f_host[0][idx];
            };
            
        } else if (device=="gpu") {
            (*f_ptr_i).f_temp_host=new double[(*f_ptr_i).gridNumberAll()];
            for (int i_f_copy=0; i_f_copy<num_field_copy; i_f_copy++) {
                cudaMalloc(&(*f_ptr_i).f[i_f_copy], num_grid*sizeof(double));
                if ((*f_ptr_i).priority()==0) {
                    cudaMalloc(&(*f_ptr_i).rhs[i_f_copy], num_grid*sizeof(double));
                    cudaMalloc(&(*f_ptr_i).lhs[i_f_copy], num_grid*sizeof(double));
                };
            };
            cudaMemcpy((*f_ptr_i).f[0], (*f_ptr_i).f_host[0], num_grid*sizeof(double),cudaMemcpyHostToDevice);            
        };
    };
    
};


// -------------------------------------------------------------
void Evolver::initRHSs() {
    
    // Loop over fields to get all functions appearing in their RHS.
    for (auto f_ptr_i : (*system_ptr).field_ptrs ) {
        int num_grid=(*f_ptr_i).gridNumberAll();
        // Number of terms appeared on this field's RHS
        int num_terms=(*f_ptr_i).rhsTerms().size();
        // Get total number of functions
        int num_funcs=0;
        for (auto rhs_term_i : (*f_ptr_i).rhsTerms()) {
            num_funcs+=rhs_term_i.f_funcs.size();
        };
        // Number of [0] explicit and [1] implicit terms
        (*f_ptr_i).rhs_ptrs_host.num_terms=new int[2];
        (*f_ptr_i).rhs_ptrs_host.num_funcs_1term=new int[num_terms];
        // Prefactors of each term
        (*f_ptr_i).rhs_ptrs_host.prefactors=new double[num_terms];
        // Numerical schemes used to calculate functions
        (*f_ptr_i).rhs_ptrs_host.schemes=new int[num_funcs];
        // Pointers to memory storing the function fields
        (*f_ptr_i).rhs_ptrs_host.f_func_ptrs=new double * [num_funcs];
        
        
        // Loop over terms on the RHS of each field.
        int i_func=0;
        int i_term=0;
        int num_terms_expl=0;
        
        // Loop over all explicit terms to add functions
        for (auto rhs_term_i : (*f_ptr_i).rhsTerms()) {
            int i_func_1term=0;
            if (rhs_term_i.scheme=="explicit") {
                num_terms_expl+=1;
                
                // Loop over operators in each term
                for (auto f_func_i : rhs_term_i.f_funcs) {
                    // Add this term to function list of this field. Each function appears only once.                    
                    // Check if function is already in the list                    
                    int toAdd=1;
                    for (int i=0; i<(*f_func_i.field_ptr).num_f_funcs; i++) {
                        FFuncItem f_func_i1=(*f_func_i.field_ptr).f_funcs_rhs[i];

                        // If this function is already in the list
                        if (f_func_i.f_operator == f_func_i1.f_operator) {
                            toAdd=0;
                        };                        
                    };
                    if (toAdd==1) {                        
                        (*f_func_i.field_ptr).addFunctoRHS(f_func_i,device,FDMScheme);
                    };
                    
                    // Add field function pointer to rhs_ptrs_host.
                    (*f_ptr_i).rhs_ptrs_host.f_func_ptrs[i_func]=(*f_func_i.field_ptr).getFFuncPtr(f_func_i.f_operator);

                    i_func+=1;
                    i_func_1term+=1;
                };
                (*f_ptr_i).rhs_ptrs_host.num_funcs_1term[i_term]=i_func_1term;
                (*f_ptr_i).rhs_ptrs_host.prefactors[i_term]=rhs_term_i.prefactor;
                i_term+=1;
            };
        };

        int num_terms_impl=0;
        for (auto rhs_term_i : (*f_ptr_i).rhsTerms()) {            
            if (rhs_term_i.scheme=="semiImplicit") {
                int i_func_1term=0;
                num_terms_impl+=1;
                // (*f_ptr_i).rhs_ptrs_host.schemes[i_term]=-1;                
                // Evaluate each operator applied on field
                for (auto f_func_i : rhs_term_i.f_funcs) {                    
                
                    // Add this term to function list of this field. Each function appears only once.                    
                    int toAdd=1;
                    // Check if function is already in the list
                    for (int i=0; i<(*f_func_i.field_ptr).num_f_funcs; i++) {
                        FFuncItem f_func_i1=(*f_func_i.field_ptr).f_funcs_rhs[i];
                        if (f_func_i.f_operator == f_func_i1.f_operator) {
                            toAdd=0;
                        };
                    };
                    if (toAdd==1) {                            
                        // (*f_ptr_i).f_funcs_rhs.push_back(f_func_i);
                        (*f_func_i.field_ptr).addFunctoRHS(f_func_i,device,FDMScheme);
                    };

                    (*f_ptr_i).rhs_ptrs_host.f_func_ptrs[i_func]=(*f_func_i.field_ptr).getFFuncPtr(f_func_i.f_operator);

                    i_func+=1;
                    i_func_1term+=1;
                };
                (*f_ptr_i).rhs_ptrs_host.num_funcs_1term[i_term]=i_func_1term;
                (*f_ptr_i).rhs_ptrs_host.prefactors[i_term]=rhs_term_i.prefactor;
                i_term+=1;
            };
        };
        
        (*f_ptr_i).rhs_ptrs_host.num_terms[0]=num_terms_expl;
        (*f_ptr_i).rhs_ptrs_host.num_terms[1]=num_terms_impl;

        // Copy these values to device.
        if (device=="gpu") {                        
            cudaMalloc(&(*f_ptr_i).rhs_ptrs_dev.num_terms,2*sizeof(int));
            cudaMalloc(&(*f_ptr_i).rhs_ptrs_dev.num_funcs_1term,num_terms*sizeof(int));
            cudaMalloc(&(*f_ptr_i).rhs_ptrs_dev.schemes,num_terms*sizeof(int));
            cudaMalloc(&(*f_ptr_i).rhs_ptrs_dev.prefactors,num_terms*sizeof(double));
            cudaMalloc(&(*f_ptr_i).rhs_ptrs_dev.f_func_ptrs,num_funcs*sizeof(double*));
            cudaMalloc(&(*f_ptr_i).f_funcs_dev,200*sizeof(double*));
            cudaMemcpy((*f_ptr_i).rhs_ptrs_dev.num_terms, (*f_ptr_i).rhs_ptrs_host.num_terms, 2*sizeof(int),cudaMemcpyHostToDevice);
            cudaMemcpy((*f_ptr_i).rhs_ptrs_dev.num_funcs_1term, (*f_ptr_i).rhs_ptrs_host.num_funcs_1term, num_terms*sizeof(int),cudaMemcpyHostToDevice);
            cudaMemcpy((*f_ptr_i).rhs_ptrs_dev.schemes, (*f_ptr_i).rhs_ptrs_host.schemes, num_terms*sizeof(int),cudaMemcpyHostToDevice);
            cudaMemcpy((*f_ptr_i).rhs_ptrs_dev.prefactors, (*f_ptr_i).rhs_ptrs_host.prefactors, num_terms*sizeof(double),cudaMemcpyHostToDevice);
            cudaMemcpy((*f_ptr_i).rhs_ptrs_dev.f_func_ptrs, (*f_ptr_i).rhs_ptrs_host.f_func_ptrs, num_funcs*sizeof(double*),cudaMemcpyHostToDevice);
            cudaMemcpy((*f_ptr_i).f_funcs_dev, (*f_ptr_i).f_funcs_host, 200*sizeof(double*),cudaMemcpyHostToDevice);
        };

    };    
};


// =============================================================

#endif

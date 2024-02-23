#!/usr/bin/perl -w

# system("nvcc -I /usr/include/cuda/ -L /usr/lib64/ --cudart=shared -lcufft active_nematic_droplet_grow.cu -o active_nematic_droplet_grow");
# system("nvcc -I /usr/include/cuda/ -L /usr/lib64/ --cudart=shared -lcufft active_nematic_droplet_test1.cu -o active_nematic_droplet_test");
system("nvcc -I /usr/include/cuda/ -L /usr/lib64/ --cudart=shared -lcufft IncompFlowTest.cu -o poissonTest");

# if (-d "data") {
#     system("rm data/*");
# } else {
#     system("mkdir data");
# }

system("./poissonTest");

#!/usr/bin/perl -w

# system("nvcc -I /usr/include/cuda/ -L /usr/lib64/ --cudart=shared -lcufft active_nematic_droplet_grow.cu -o active_nematic_droplet_grow");
system("nvcc -I /usr/include/cuda/ -L /usr/lib64/ --cudart=shared -lcufft active_nematic_droplet_test.cu -o active_nematic_droplet_test");

if (-d "data") {
    system("rm data/*");
} else {
    system("mkdir data");
}

system("./active_nematic_droplet_test");

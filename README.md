# dimice-heat-dealii
# Build
$ cd dimice-heat-dealii  
$ mkdir build  
$ cmake ..  
$ make
# Examples
## Run with default parameters
$ cd dimice-heat-dealii  
$ mkdir run  
$ cd run  
$ ../build/heat
## Run with a parameter input file
$ ../build/heat ../tests/poisson_problem_01.prm

# dimice-heat-dealii
Compile this by running

	git clone alexanderzimmerman/dimice-heat-dealii
	mkdir build  
	cmake ..  
	make
	
# Examples
These examples assumes that some environment variables have been set, for example

	export DIMICE_HEAT_ROOT=~/dimice-heat-dealii
	export DIMICE_HEAT_EXE=$DIMICE_HEAT_ROOT/build/heat

Run with default parameters by running

	mkdir run  
	cd run  
	$DIMICE_HEAT_EXE
	
Run with a parameter input file by running, for example

	$DIMICE_HEAT_EXE $DIMICE_HEAT_ROOT/inputs/sphere.prm

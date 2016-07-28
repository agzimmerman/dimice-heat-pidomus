# dimice-heat-dealii
First pi-DoMUS must be installed.

Set an environment variable pointing to pi-DoMUS, e.g.

	export PIDOMUS_ROOT=/usr/local/pidomus

Compile dimice-heat-dealii by running

	git clone git@github.com:alexanderzimmerman/dimice-heat-dealii.git
	mkdir build
	cd build
	cmake -DPIDOMUS_DIR=$PIDOMUS_DIR -DCMAKE_INSTALL_PREFIX=/usr/local/dimice-heat .. 
	sudo make install
	
# Examples
These examples assumes that some environment variables have been set, for example

	export DIMICE_HEAT_ROOT=~/dimice-heat-dealii
	alias dimice-heat=$DIMICE_HEAT_ROOT/build/heat

Run with default parameters by running

	mkdir run  
	cd run  
	dimice-heat
	
Run with a parameter input file by running, for example

	dimice-heat $DIMICE_HEAT_ROOT/inputs/sphere-cylinder.prm

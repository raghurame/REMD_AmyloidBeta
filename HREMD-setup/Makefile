all:
        # ~~~~~~~~~~~~~~~~~
        # AVAILABLE OPTIONS
        # ~~~~~~~~~~~~~~~~~
	# 1. Type 'make compile' to compile the code
	# 2. Type 'make run' to run the code.
compile:
	gcc -o modifyHamiltonian modifyHamiltonian.c -lm -Wall
run:
	./modifyHamiltonian ffbonded.itp ffnonbonded.itp topol.top ffbonded_modified.itp ffnonbonded_modified.itp topol_modified.top

all:
        ###################
        # AVAILABLE OPTIONS
        ###################
	# 1. Type 'make compile' to compile the code
	# 2. Type 'make run' to run the code.
compile:
	gcc -o modifyHamiltonian modifyHamiltonian.c -lm -Wall -fstack-protector
run:
	./modifyHamiltonian ffbonded.itp ffnonbonded.itp topol.top ffbonded_0.95.itp ffnonbonded_0.95.itp topol_0.95.top 0.95 hotResidues.input hotInteractions.input
	./modifyHamiltonian ffbonded.itp ffnonbonded.itp topol.top ffbonded_0.9.itp ffnonbonded_0.9.itp topol_0.9.top 0.9 hotResidues.input hotInteractions.input
	./modifyHamiltonian ffbonded.itp ffnonbonded.itp topol.top ffbonded_0.85.itp ffnonbonded_0.85.itp topol_0.85.top 0.85 hotResidues.input hotInteractions.input
	./modifyHamiltonian ffbonded.itp ffnonbonded.itp topol.top ffbonded_0.8.itp ffnonbonded_0.8.itp topol_0.8.top 0.8 hotResidues.input hotInteractions.input
	./modifyHamiltonian ffbonded.itp ffnonbonded.itp topol.top ffbonded_0.75.itp ffnonbonded_0.75.itp topol_0.75.top 0.75 hotResidues.input hotInteractions.input
	./modifyHamiltonian ffbonded.itp ffnonbonded.itp topol.top ffbonded_0.7.itp ffnonbonded_0.7.itp topol_0.7.top 0.7 hotResidues.input hotInteractions.input
	./modifyHamiltonian ffbonded.itp ffnonbonded.itp topol.top ffbonded_0.65.itp ffnonbonded_0.65.itp topol_0.65.top 0.65 hotResidues.input hotInteractions.input
	./modifyHamiltonian ffbonded.itp ffnonbonded.itp topol.top ffbonded_0.6.itp ffnonbonded_0.6.itp topol_0.6.top 0.6 hotResidues.input hotInteractions.input
	./modifyHamiltonian ffbonded.itp ffnonbonded.itp topol.top ffbonded_0.55.itp ffnonbonded_0.55.itp topol_0.55.top 0.55 hotResidues.input hotInteractions.input
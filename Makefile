all:
	# ~~~~~~~~~~~~~~~~~
	# AVAILABLE OPTIONS
	# ~~~~~~~~~~~~~~~~~
	# 1. Type "make pedist" to compile distribution c file.
pedist:
	gcc -o PEdistribution PEdistribution.c -lm -Wall

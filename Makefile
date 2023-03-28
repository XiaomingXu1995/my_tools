app=simulate-mutation-sequence cmp-dist-RabbitKSSD-bindash
CXXFLAGS=-O3 -g -Wall
all: ${app}
simulate-mutation-sequence: src/simulate_mutation_sequence.cpp
	g++ ${CXXFLAGS} src/simulate_mutation_sequence.cpp -o simulate-mutation-sequence
cmp-dist-RabbitKSSD-bindash: src/cmp_dist_RabbitKSSD_binDash.cpp
	g++ ${CXXFLAGS} -fopenmp src/cmp_dist_RabbitKSSD_binDash.cpp -o cmp-dist-RabbitKSSD-bindash

clean:
	rm ${app}
	

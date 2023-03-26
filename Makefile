app=simulate-mutation-sequence
CXXFLAGS=-O3 -g -Wall
all: ${app}
simulate-mutation-sequence: src/simulate_mutation_sequence.cpp
	g++ ${CXXFLAGS} src/simulate_mutation_sequence.cpp -o simulate-mutation-sequence

clean:
	rm ${app}
	

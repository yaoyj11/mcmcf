CC=g++
FLAGS=-O3 -std=c++11
LIB = -lemon

all: mcmcf

mcmcf:flow_solution.o main.o fractional_packing.o
	$(CC) $(FLAGS) $(LIB) -o mcmcf main.o flow_solution.o fractional_packing.o

flow_solution.o: flow_solution.cc demand.h flow_solution.h
	$(CC) $(FLAGS) -c flow_solution.cc

fractional_packing.o: fractional_packing.cc fractional_packing.h flow_solution.h demand.h
	$(CC) $(FLAGS) -c fractional_packing.cc

main.o: main.cc fractional_packing.h
	$(CC) $(FLAGS) -c main.cc

clean:
	rm *.o
	rm mcmcf

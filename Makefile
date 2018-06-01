CC=g++
FLAGS=-std=c++11 -O3
LIB = -lemon

all: mcmcf test test-or
	./test
	./test-or

mcmcf:flow_solution.o main.o fractional_packing.o
	$(CC) $(FLAGS) $(LIB) -o mcmcf main.o flow_solution.o fractional_packing.o

test:flow_solution.o test.o fractional_packing.o
	$(CC) $(FLAGS) $(LIB) -o test test.o flow_solution.o fractional_packing.o

test-or: test-or.o
	$(CC) $(FLAGS) $(LIB) -lortools -o test-or test-or.o

flow_solution.o: flow_solution.cc demand.h flow_solution.h
	$(CC) $(FLAGS) -c flow_solution.cc

fractional_packing.o: fractional_packing.cc fractional_packing.h flow_solution.h demand.h
	$(CC) $(FLAGS) -c fractional_packing.cc

main.o: main.cc fractional_packing.h
	$(CC) $(FLAGS) -c main.cc

test.o: test.cc fractional_packing.h
	$(CC) $(FLAGS) -c test.cc

test-or.o: test-or.cc
	$(CC) $(FLAGS) -c test-or.cc
clean:
	rm *.o
	rm mcmcf
	rm test

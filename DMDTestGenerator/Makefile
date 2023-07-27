CC=g++ -DOzgurTest
# flags -DKutzTest

all: solver.x
	rm -rf *.o

solver.x: dmdInputDataGenerator.o
	$(CC) dmdInputDataGenerator.o -o solver.x

dmdInputDataGenerator.o: dmdInputDataGenerator.cpp
	$(CC) -c --std=c++11 dmdInputDataGenerator.cpp

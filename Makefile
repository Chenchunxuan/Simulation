SimObjects = main.o 
TestObjects = test_funcs.o


sim: $(SimObjects)
	g++ $(SimObjects) -o sim

test_funcs: $(TestObjects)
	g++ $(TestObjects) -o test_funcs


main.o: main.cpp main.hpp matrix.hpp
	g++ -c main.cpp 

test_funcs.o: test_funcs.cpp 
	g++ -c test_funcs.cpp

clean:
	rm *.o sim test_funcs
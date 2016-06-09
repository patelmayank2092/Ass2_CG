
CDIR = myColsamm/Source/
CFLAG = -c -O3 -std=c++11 
FFLAG = -O3 -std=c++11 -I $(CDIR)



all:	main initialisation 

#solve:	initialisation.o matrix_generator main.o
#	g++ -O3 -std=c++11 initialisation.o matrix_generator.o main.o -o solve

initialisation:   initialisation.h initialisation.cpp
	g++ $(CFLAG) initialisation.cpp -o initialisation.o

main:   type.h  Vertex.h initialisation.h matrix_generator.h  $(CDIR)*.h main.cpp
	g++ $(FFLAG) initialisation.cpp matrix_generator.cpp main.cpp -o result

#matrix_generator: matrix_generator.h matrix_generator.cpp 
#	g++ -c -O3 -std=c++11 matrix_generator.cpp -o matrix_generator.o

clean:
	rm -rf *.o
	rm result

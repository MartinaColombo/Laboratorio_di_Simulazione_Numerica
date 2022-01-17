MPI_COMPILE_FLAG = $(shell mpicc --showme:compile)
MPI_LINK_FLAGS = $(shell mpicc --showme:link)
CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o funzioni.o random.o
	$(CC) $(MPI_COMPILE_FLAGS) funzioni.o random.o main.o -o $(MPI_LINK_FLAGS) -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
clean :
	rm *.o main.exe L2_rank1.dat L2_rank2.dat L2_rank3.dat L2_rank4.dat aveL2_rank1.dat aveL2_rank2.dat aveL2_rank3.dat aveL2_rank4.dat Best_crom_rank1.dat Best_crom_rank2.dat Best_crom_rank3.dat Best_crom_rank4.dat

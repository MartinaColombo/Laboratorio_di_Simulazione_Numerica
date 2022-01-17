CC = g++
CFLAGS = -Wall -O3 --std=c++11

MolDyn_NVE.exe : MolDyn_NVE.o 
	$(CC)  MolDyn_NVE.o -o MolDyn_NVE.exe
MolDyn_NVE.o : MolDyn_NVE.cpp MolDyn_NVE.h
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
clean :
	rm *.o MolDyn_NVE.exe output.ekin.0 output.epot.0 output.etot.0 output.gave.0 output.gofr.0 output.temp.0 output_ekin.dat output_epot.dat output_etot.dat output_temp.dat output.press.0 output_press.dat

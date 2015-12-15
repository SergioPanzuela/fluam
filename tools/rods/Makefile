CC=g++

#DEBUG=-g -DDEBUG

INCLUDES = 



CFLAGS= -Ofast -ffast-math $(DEBUG) -funroll-loops  -march=native -std=c++11 -pipe -fomit-frame-pointer -funroll-loops



EXECUTABLE=mesh

all: mesh rods

mesh: mesh.o 
	$(CC) -o mesh mesh.o
rods: rods.o 
	$(CC) -o rods rods.o
.cpp.o: 
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@



clean: 
	rm mesh.o rods.o mesh rods
redo: clean all

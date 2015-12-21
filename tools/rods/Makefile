CC=clang++-3.6

#DEBUG=-g -DDEBUG

OBJS=mesh.o rods.o 

CFLAGS= -Ofast -ffast-math $(DEBUG) -funroll-loops  -march=native -std=c++11 -pipe -fomit-frame-pointer -funroll-loops
CXXFLAGS=$(CFLAGS)

EXECUTABLE=mesh

all: mesh rods

clean: 
	rm -f $(OBJS) mesh rods
redo: clean all

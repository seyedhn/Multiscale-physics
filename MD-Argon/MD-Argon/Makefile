BASE = /Users/Seyed/Documents/Startup/C++/MD-Argon
CC = g++
CFLAGS = -c

INCLS = -I/Users/Seyed/Documents/Startup/C++/eigen/ -I/Users/Seyed/Documents/Startup/C++/gnuplot-iostream/ -I/Users/Seyed/Documents/Startup/C++/boost/


EFILE = pgm

all: pgm 

pgm: main.o
	$(CC) $(INCLS) main.o -o $(EFILE) 

main.o: main.cpp
	$(CC) $(INCLS) main.cpp 

clean:
	rm -rf *o pgm
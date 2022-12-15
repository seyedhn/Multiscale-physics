BASE = /Users/Seyed/Documents/Startup/C++/MC-Ising
CC = c++
CFLAGS = -c  -Wall
INCLS = -I/usr/local/include -I/Users/Seyed/Documents/Startup/C++/eigen/ -I/Users/Seyed/Documents/Startup/C++/gnuplot-iostream/ -lgmp -lgmpxx

EFILE = pgm

all: pgm 

pgm: main.o MC.o
	$(CC) main.o MC.o -o $(EFILE)  $(INCLS)

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp -o main.o   $(INCLS)	


MC.o: MC.cpp
	$(CC) $(CFLAGS) MC.cpp   $(INCLS)

clean:
	rm -rf *o pgm
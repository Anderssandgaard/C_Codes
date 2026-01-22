CFLAGS +=  -Wall -std=gnu1x $$(gsl-config --cflags) -fopenmp
LDLIBS += -lm -lgomp $$(gsl-config --libs)

.PHONEY: all
all: out.txt 

out.txt: main Makefile
	./$< 1>$@

main: main.o Monte_Carlo.o matrixandvectors.o
Monte_Carlo: Monte_Carlo.o matrixandvectors.o


main.o Monte_Carlo.o matrixandvectors.o: matrixandvectors.h

matrixandvectors.o: matrixandvectors.c

clean: ; $(RM)\
		main *.o \
		Monte_Carlo *.o\




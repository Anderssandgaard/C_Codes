CFLAGS +=  -Wall -std=gnu1x 
LDLIBS += -lm 

.PHONEY: all
all: out.txt 

out.txt: main Makefile
	./$< 1>$@ 

main: main.o Integration.o matrixandvectors.o
Integration: Integration.o matrixandvectors.o

main.o Integration.o matrixandvectors.o: matrixandvectors.h

matrixandvectors.o: matrixandvectors.c

clean: ; $(RM)\
		main *.o \
		Integration *.o\


CFLAGS = -Wall -std=gnu1x
LDLIBS = -lm


out.txt: main
	./$< > $@ 

main: main.o QR_decomp.o matrixandvectors.o
QR_decomp: QR_decomp.o matrixandvectors.o

main.o QR_decomp.o matrixandvectors.o: matrixandvectors.h
matrixandvectors.o: matrixandvectors.c

clean: ; $(RM) \
		main *.o \
		QR_decomp *.o \
		matrixandvectors.o
		


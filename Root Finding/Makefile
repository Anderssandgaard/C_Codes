CFLAGS +=  -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += -lm $$(gsl-config --libs)

.PHONEY: all
all: out.txt out.gsl.Rosenbrock.txt  out.gsl.Himmelblau.txt out.gsl.Sys.txt

out.txt: main Makefile
	./$< >$@ 

out.gsl.Sys.txt: main_gsl_Sys Makefile
	./$< >$@ 

out.gsl.Rosenbrock.txt: main_gsl_Rosenbrock Makefile
	./$< >$@ 

out.gsl.Himmelblau.txt: main_gsl_Himmelblau Makefile
	./$< >$@ 

main: main.o Root_Finding.o matrixandvectors.o
Root_Finding: Root_Finding.o matrixandvectors.o
main_gsl_Rosenbrock: main_gsl_Rosenbrock.o 
main_gsl_Himmelblau: main_gsl_Himmelblau.o 
main_gsl_Sys: main_gsl_Sys.o 

main.o Root_Finding.o matrixandvectors.o: matrixandvectors.h

matrixandvectors.o: matrixandvectors.c

clean: ; $(RM)\
		main *.o \
		main_gsl_Rosenbrock *.o \
		main_gsl_Himmelblau *.o \
		main_gsl_Sys *.o \
		Root_Finding *.o\



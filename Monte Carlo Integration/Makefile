CFLAGS +=  -Wall -std=gnu1x 
LDLIBS += -lm 

.PHONEY: all
all: out.txt Monte_carlo_error.svg

Monte_carlo_error.svg: plot.gpi datalog
	gnuplot $<

out.txt: main Makefile
	./$< 1>$@ 2>datalog

main: main.o Monte_Carlo.o matrixandvectors.o
Monte_Carlo: Monte_Carlo.o matrixandvectors.o


main.o Monte_Carlo.o matrixandvectors.o: matrixandvectors.h

matrixandvectors.o: matrixandvectors.c

clean: ; $(RM)\
		main *.o \
		Monte_Carlo *.o\

plot.gpi: Makefile
	echo 'set term svg size 800,600 background rgb "white"' > $@
	echo 'set out  "Monte_carlo_error.svg"'>> $@
	echo 'set xlabel "N(1)"' >> $@
	echo 'set ylabel "Error(1)"' >> $@
	echo 'set key right bottom font ",10"' >> $@
	echo 'set title  "Error of Plain Monte Carlo Integration on Func. x*x+y*y+z*z"' >> $@
	echo 'f(x)=b+(a/sqrt(x))' >>$@
	echo 'b=0.1; a=10' >>$@
	echo 'fit f(x) "datalog" via b,a' >>$@
	echo 'F = sprintf("fit: error=%4.1f+(%3.0f/sqrt(N))",b,a)' >>$@
	echo 'plot "datalog" using 1:2 with points  title "Data of Error vs. N",f(x) with lines title F \' >> $@
	




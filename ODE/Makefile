CFLAGS +=  -Wall -std=gnu1x 
LDLIBS += -lm 

.PHONEY: all
all: out.txt ODE.svg

ODE.svg: plot.gpi datalog
	gnuplot $<

out.txt: main Makefile
	./$< 1>$@ 2>datalog

main: main.o ODE.o matrixandvectors.o
ODE: ODE.o matrixandvectors.o


main.o ODE.o matrixandvectors.o: matrixandvectors.h

matrixandvectors.o: matrixandvectors.c

clean: ; $(RM)\
		main *.o \
		ODE *.o\

plot.gpi: Makefile
	echo 'set term svg size 800,600 background rgb "white"' > $@
	echo 'set out  "ODE.svg"'>> $@
	echo 'set xlabel "x"' >> $@
	echo 'set ylabel "y"' >> $@
	echo 'set key right bottom font ",10"' >> $@
	echo 'set title  "ODE for sin(x) and 1/(1+exp(-x))"' >> $@
	echo 'plot "datalog" index 0 using 1:2 with lines  title "ODE RK23 Harmonic Diff. Eqn." \' >> $@
	echo ',"datalog" index 0 using 1:3 with lines  title "Theoretical Sine Function"\' >> $@
	echo ',"datalog" index 1 using 1:2 with lines  title "ODE RK23 Logistic Diff. Eqn." \' >> $@
	echo ',"datalog" index 1 using 1:3 with lines title "Theoretical Log. Func."' >> $@




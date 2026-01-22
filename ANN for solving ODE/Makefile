CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

Neural_Network_ODE.svg: plot.gpi out.txt Makefile
	gnuplot $<

out.txt:main
	./$< 1> $@ 2>log

main:main.o neurons.o

clean:
	rm -f main *.o neurons.o

plot.gpi: Makefile
	echo 'set term svg size 800,600 background rgb "white"' > $@
	echo 'set out  "Neural_Network_ODE.svg"'>> $@
	echo 'set xlabel "x"' >> $@
	echo 'set ylabel "y"' >> $@
	echo 'set key right bottom font ",10"' >> $@
	echo 'set title  "Neural Network for solving ODE "' >> $@
	echo 'plot "out.txt" index 0 with lines title "Neural Logistic Function" \' >> $@
	echo ', 1/(1+exp(-x)) with lines title "Analytic Logistic Function"\' >> $@
	echo ', exp(-x*x*0.5) with lines title "Analytic Gaussian Function"\' >> $@
	echo ',"out.txt" index 1 with lines title "Neural Gaussian Function"' >> $@

CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

Neural_Network_Interpolation.svg: plot.gpi out.txt Makefile
	gnuplot $<

out.txt:main
	./$< 1> $@ 2>log

main:main.o neurons.o

clean:
	rm -f main *.o neurons *.o

plot.gpi: Makefile
	echo 'set term svg size 800,600 background rgb "white"' > $@
	echo 'set out  "Neural_Network_Interpolation.svg"'>> $@
	echo 'set xlabel "x"' >> $@
	echo 'set ylabel "y"' >> $@
	echo 'set key right bottom font ",10"' >> $@
	echo 'set title  "Neural Network Interpolation of cos(5*x-1)*exp(-x*x) "' >> $@
	echo 'plot "out.txt" index 0 with points title "Data","out.txt" index 1 with lines title "Neural Interpolation" ' >> $@
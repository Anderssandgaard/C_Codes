CFLAGS = -Wall -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

.PHONEY: default
default: plot.pdf

plot.pdf: plot.gpi out.txt 
	gnuplot $<

out.txt: main
	./$< $$(wc -l points.data) < points.data 1>$@ 
main: main.o lininterp.o quadinterp.o cubicinterp.o

clean: ; $(RM)\
		main *.o \
		 lininterp *.o \
		 quadinterp *.o \
		 cubicinterp *.o \
		


plot.gpi: Makefile
	echo 'set term pdf' > $@
	echo 'set out  "plot.pdf"'>> $@
	echo 'set xlabel "x"' >> $@
	echo 'set ylabel "y"' >> $@
	echo 'set grid' >> $@
	echo 'set key right bottom' >> $@
	echo 'set title  "Spline of Random Points"' >> $@
	echo 'plot "out.txt" using 1:2 with lines title "Linear Spline"\' >> $@
	echo ',"out.txt" using 1:4 with lines title "Quad Spline"\'>>  $@
	echo ',"out.txt" using 1:12 with lines title "Cubic Spline"\'>>  $@
	echo ',"out.txt" using 1:7 with lines title "GSL Lin. Spline"\'  >> $@
	echo ',"out.txt" using 1:9 with lines title "GSL Cubic. Spline"'  >> $@
	echo 'set title  "Ïntegration Spline of Random Points"' >> $@
	echo 'plot "out.txt" using 1:3 with lines title "Linear Spline integ"\' >> $@
	echo ',"out.txt" using 1:5 with lines title "Quad Spline integ"\' >> $@
	echo ',"out.txt" using 1:13 with lines title "Cubic Spline integ"\' >> $@
	echo ',"out.txt" using 1:8 with lines title "GSL Lin. integ Spline"\' >> $@
	echo ',"out.txt" using 1:10 with lines title "GSL Cubic. integ Spline"' >> $@
	echo 'set title  "Derivative Spline of Random Points"' >> $@
	echo 'plot "out.txt" using 1:6 with lines title "Quad Spline deriv"\' >> $@
	echo ',"out.txt" using 1:14 with lines title "Cubic Spline deriv"\' >> $@
	echo ',"out.txt" using 1:11 with lines title "GSL Cubic. deriv Spline"' >> $@

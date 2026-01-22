CFLAGS +=  -Wall -std=gnu1x
LDLIBS += -lm


.PHONEY: all
all: Least-Squares.pdf out.txt out.svd.txt

Least-Squares.pdf: plot.gpi data.txt out.txt out.svd.txt
	gnuplot $<

out.txt: main Makefile
	./$< $$(wc -l data.txt) < data.txt 1>$@ 

out.svd.txt: main_svd Makefile
	./$< $$(wc -l data.txt) < data.txt 1>$@ 

main: main.o Least-Squares.o matrixandvectors.o
main_svd: main_svd.o Least-Squares.o matrixandvectors.o
Least-Squares: Least-Squares.o matrixandvectors.o

main.o main_svd.o Least-Squares.o matrixandvectors.o: matrixandvectors.h
matrixandvectors.o: matrixandvectors.c

clean: ; $(RM)\
		main *.o \
		main_svd *.o \
		Least-Squares.o \
		matrixandvectors.o \

plot.gpi: Makefile
	echo 'set term pdf' > $@
	echo 'set out  "Least-Squares.pdf"'>> $@
	echo 'set xlabel "x"' >> $@
	echo 'set ylabel "y"' >> $@
	echo 'set key right bottom font ",10"' >> $@
	echo 'set title  "Least Squares Fit to c1*Log(x)+c2+c3*x for Data"' >> $@
	echo 'plot "data.txt" with errorbars title "Data"\' >> $@
	echo ',"out.txt" index 1 using 1:2 with lines title "Least Squares Fit"\' >> $@
	echo ',"out.txt" index 1 using 1:3 with lines title "Least Squares Fit plus log(x)-delta" \' >> $@
	echo ',"out.txt" index 1 using 1:4 with lines title "Least Squares Fit minus log(x)-delta"\' >> $@
	echo ',"out.txt" index 2 using 1:3 with lines title "Least Squares Fit plus const-delta" \' >> $@
	echo ',"out.txt" index 2 using 1:4 with lines title "Least Squares Fit minus const-delta"\' >> $@
	echo ',"out.txt" index 3 using 1:3 with lines title "Least Squares Fit plus x-delta" \' >> $@
	echo ',"out.txt" index 3 using 1:4 with lines title "Least Squares Fit minus x-delta"' >> $@
	echo 'set title  "Least Squares Fit to c1*Log(x)+c2+c3*x for Data Using SVD"' >> $@
	echo 'plot "data.txt" with errorbars title "Data"\' >> $@
	echo ',"out.svd.txt" index 1 using 1:2 with lines title "Least Squares Fit"\' >> $@
	echo ',"out.svd.txt" index 1 using 1:3 with lines title "Least Squares Fit plus log(x)-delta" \' >> $@
	echo ',"out.svd.txt" index 1 using 1:4 with lines title "Least Squares Fit minus log(x)-delta"\' >> $@
	echo ',"out.svd.txt" index 2 using 1:3 with lines title "Least Squares Fit plus const-delta" \' >> $@
	echo ',"out.svd.txt" index 2 using 1:4 with lines title "Least Squares Fit minus const-delta"\' >> $@
	echo ',"out.svd.txt" index 3 using 1:3 with lines title "Least Squares Fit plus x-delta" \' >> $@
	echo ',"out.svd.txt" index 3 using 1:4 with lines title "Least Squares Fit minus x-delta"' >> $@



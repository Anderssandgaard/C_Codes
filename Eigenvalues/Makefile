CFLAGS +=  -std=gnu1x
CFLAGS += -march=native -Ofast
LDLIBS += -lm


.PHONEY: all
all: times.pdf out.txt

times.pdf: times.gpi out.times.rowbyrow.txt out.times.classic.txt out.times.eigbyeig.txt out.times.lowesteig.txt
	gnuplot $<

out.txt: main Makefile
	./$< 5 > $@ 

out.times.rowbyrow.txt: main_rbr Makefile
	cat /dev/null > $@
	for n in `seq 150 2 200`; do time --format "$$n %U" --append --output $@ ./main_rbr $$n 2>/dev/null; done

out.times.eigbyeig.txt: main_eigbyeig Makefile
	cat /dev/null > $@
	for n in `seq 50 2 100`; do time --format "$$n %U" --append --output $@ ./main_eigbyeig $$n 2>/dev/null; done

out.times.lowesteig.txt: main_lowesteig Makefile
	cat /dev/null > $@
	for n in `seq 150 2 200`; do time --format "$$n %U" --append --output $@ ./main_lowesteig $$n 2>/dev/null; done

out.times.classic.txt: main_classic Makefile
	cat /dev/null > $@
	for n in `seq 150 2 200`; do time --format "$$n %U" --append --output $@ ./main_classic $$n 2>/dev/null; done	



main: main.o Matrix_Diagonalisation.o matrixandvectors.o
main_rbr: main_rbr.o Matrix_Diagonalisation.o matrixandvectors.o
main_eigbyeig: main_eigbyeig.o Matrix_Diagonalisation.o matrixandvectors.o
main_lowesteig: main_lowesteig.o Matrix_Diagonalisation.o matrixandvectors.o
main_classic: main_classic.o Matrix_Diagonalisation.o matrixandvectors.o
Matrix_Diagonalisation: Matrix_Diagonalisation.o matrixandvectors.o


main.o main_classic main_rbr main_eigbyeig main_lowesteig Matrix_Diagonalisation.o matrixandvectors.o: matrixandvectors.h
matrixandvectors.o: matrixandvectors.c

clean: ; $(RM)\
		main *.o \
		main_rbr *.o \
		main_eigbyeig *.o \
		main_lowesteig *.o \
		main_classic *.o \
		Matrix_Diagonalisation *.o \
		matrixandvectors.o \
		times.gpi

times.gpi: Makefile
	echo 'set term pdf' > $@
	echo 'set out  "times.pdf"'>> $@
	echo 'set xlabel "Matrix Size"' >> $@
	echo 'set ylabel "Diagonalisation Time (s)"' >> $@
	echo 'set key left top font ",10"' >> $@
	echo 'set title  "Diagonalisation Time vs. Matrix Dimension"' >> $@
	echo 'f(x)=b+(x/a)**3' >>$@
	echo 'b=0.1; a=10' >>$@
	echo 'fit f(x) "out.times.rowbyrow.txt" via b,a' >>$@
	echo 'F = sprintf("fit: %4.1f+(n/%3.0f)^3",b,a)' >>$@
	echo 'ff(x)=c+(x/d)**3' >>$@
	echo 'c=0.1; d=10' >>$@
	echo 'fit ff(x) "out.times.classic.txt" via c,d' >>$@
	echo 'FF = sprintf("fit: %4.1f+(n/%3.0f)^3",c,d)' >>$@
	echo 'fff(x)=e+(x/f)**3' >>$@
	echo 'e=3; f=3' >>$@
	echo 'fit fff(x) "out.times.eigbyeig.txt" via e,f' >>$@
	echo 'FFF = sprintf("fit: %4.1f+(n/%3.0f)^3",e,f)' >>$@
	echo 'ffff(x)=p1+(p2/a)**3' >>$@
	echo 'p1=0.1; p2=10' >>$@
	echo 'fit ffff(x) "out.times.lowesteig.txt" via p1,p2' >>$@
	echo 'FFFF = sprintf("fit: %4.1f+(n/%3.0f)^3",p1,p2)' >>$@
	echo 'plot "out.times.rowbyrow.txt" title "Time Measurements Row-by-Row",f(x) title F \' >> $@
	echo ', "out.times.classic.txt" title "Time Measurements Classic", ff(x) title FF\' >> $@
	echo ', "out.times.lowesteig.txt" title "Time Measurements Lowest Eig", ffff(x) title FFFF' >> $@
	echo 'plot "out.times.eigbyeig.txt" title "Time Measurements Eig by Eig",fff(x) title FFF' >>$@



CFLAGS	=	-O -W -Wall -I/usr/local/include
LDFLAGS	=	-L/usr/local/lib -lmatrix -lm
PROGS	=	heat2d-i-en heat2d-i-en3 ss
OBJS	=	lu.o bandlu.o call_gnuplot.o

default: $(PROGS)

heat2d-i-en: heat2d-i-en-mk.c lu.o call_gnuplot.o
	$(CC) $(CFLAGS) -o $@ heat2d-i-en-mk.c lu.o call_gnuplot.o $(LDFLAGS)

heat2d-i-en3: heat2d-i-en3-mk.c bandlu.o call_gnuplot.o
	$(CC) $(CFLAGS) -o $@ heat2d-i-en3-mk.c bandlu.o call_gnuplot.o $(LDFLAGS)

ss: ss-mk.c bandlu.o call_gnuplot.o
	$(CC) $(CFLAGS) -o $@ ss-mk.c bandlu.o call_gnuplot.o $(LDFLAGS)

bandlu.o: bandlu.c bandlu.c
	$(CC) $(CFLAGS) -c bandlu.c

lu.o: lu.c lu.c
	$(CC) $(CFLAGS) -c lu.c

clean:
	rm -f $(PROGS) $(OBJS)

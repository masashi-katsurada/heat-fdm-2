SRCS	=	heat2d-disk-e.c heat2d-disk-i.c heat2d-disk-adi.c \
		call_gnuplot.c call_gnuplot.h ptrilu.c ptrilu.h \
		find_Bessel_roots.c
PROGS	=	bisection heat2d-disk-e heat2d-disk-i heat2d-disk-adi \
		find_Bessel_roots
CC	=	gcc
CFLAGS	=	-W -Wall -O -I/usr/local/include
LDFLAGS	=	-L/usr/local/lib -L/opt/X11/lib -lglscd -lX11 -lmatrix -lm
#LDFLAGS	=	-L/usr/local/lib -L/usr/local/X11R6.3/lib -lglscd -lX11 -lsocket -lmatrix -lm
#LDFLAGS	=	-L/usr/local/lib -L/usr/X11R6/lib -lglscd -lX11 -lmatrix -lm
OBJS1	=	heat2d-disk-e.o call_gnuplot.o ptrilu.o
OBJS2	=	heat2d-disk-i.o call_gnuplot.o ptrilu.o
OBJS3	=	heat2d-disk-adi.o call_gnuplot.o ptrilu.o trid-lu.o
OBJS	=	$(OBJS1) $(OBJS2) $(OBJS3)

.c:
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

default: $(PROGS)

heat2d-disk-e: $(OBJS1)
	$(CC) $(CFLAGS) -o $@ $(OBJS1) $(LDFLAGS)

heat2d-disk-i: $(OBJS2)
	$(CC) $(CFLAGS) -o $@ $(OBJS2) $(LDFLAGS)

heat2d-disk-adi: $(OBJS3)
	$(CC) $(CFLAGS) -o $@ $(OBJS3) $(LDFLAGS)

find_Bessel_roots: find_Bessel_roots.c
	$(CC) $(CFLAGS) -o $@ find_Bessel_roots.c -lm

expand:
	expand heat2d-disk-e.c > foo
	mv foo heat2d-disk-e.c
	expand heat2d-disk-i.c > foo
	mv foo heat2d-disk-i.c
	expand heat2d-disk-adi.c > foo
	mv foo heat2d-disk-adi.c

copy:
	scp -p $(SRCS) mathweb:web/program/fdm

clean vc:
	rm -f $(PROGS) $(OBJS)

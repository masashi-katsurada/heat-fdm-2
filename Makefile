HOST	=	nalab.mind.meiji.ac.jp
DIR	=	labo/text
DEST	=	$(HOST):Sites/$(DIR)
PUSH	=	bin/pushtosakura
FILES	=	heat-fdm-2.dvi heat-fdm-2.pdf \
		heat-fdm-2-prog.tar.gz
PROGSRC	=	prog/heat2d-disk-e.c prog/heat2d-disk-i.c \
		prog/heat2d-disk-adi.c prog/ptrilu.c Bessel/draw-Bessel.c \
		prog/find_Bessel_roots.c

default: all
all: $(PROGS) $(FILES)

heat-fdm-2.dvi: heat-fdm-2.tex ../reference/reference.tex $(PROGSRC)
	platex heat-fdm-2.tex
	-mendex heat-fdm-2
	pbibtex heat-fdm-2
	platex heat-fdm-2.tex
	platex heat-fdm-2.tex

heat-fdm-2.pdf: heat-fdm-2.dvi
	dvipdfmx -d 5 -O 2 heat-fdm-2.dvi

copy: $(FILES)
	scp -p $(FILES) $(DEST)
	ssh $(HOST) $(PUSH) $(DIR)

heat-fdm-2-prog.tar.gz:
	rm -rf heat-fdm-2-prog
	mkdir heat-fdm-2-prog
	cp -p prog/* heat-fdm-2-prog
	gtar cfz $@ heat-fdm-2-prog

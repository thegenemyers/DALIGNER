CFLAGS = -O3 -Wall -Wextra -fno-strict-aliasing

ALL = daligner HPCdaligner HPCmapper LAsort LAmerge LAsplit LAcat LAshow LAcheck

all: $(ALL)

LIB_ALIGN_HEADERS=align.h DB.h QV.h

libalign.a: align.c DB.c QV.c ${LIB_ALIGN_HEADERS}
	$(CC) $(CFLAGS) -c -o align.o align.c
	$(CC) $(CFLAGS) -c -o DB.o DB.c
	$(CC) $(CFLAGS) -c -o QV.o QV.c
	ar r libalign.a align.o DB.o QV.o
	rm -f align.o DB.o QV.o

libalign.so: align.c DB.c QV.c ${LIB_ALIGN_HEADERS}
	$(CC) $(CFLAGS) -c -fpic -o align.so align.c
	$(CC) $(CFLAGS) -c -fpic -o DB.so DB.c
	$(CC) $(CFLAGS) -c -fpic -o QV.so QV.c
	$(CC) -shared -o libalign.so align.so DB.so QV.so
	rm -f align.so DB.so QV.so

daligner: daligner.c filter.c filter.h ${LIB_ALIGN_HEADERS} libalign.a
	$(CC) $(CFLAGS) -o daligner daligner.c filter.c libalign.a -lpthread -lm 

HPCdaligner: HPCdaligner.c DB.h QV.h libalign.a
	$(CC) $(CFLAGS) -o HPCdaligner HPCdaligner.c libalign.a -lm

HPCmapper: HPCmapper.c DB.h QV.h libalign.a
	$(CC) $(CFLAGS) -o HPCmapper HPCmapper.c libalign.a -lm

LAsort: LAsort.c align.h DB.h QV.h
	$(CC) $(CFLAGS) -o LAsort LAsort.c libalign.a -lm

LAmerge: LAmerge.c align.h DB.h QV.h
	$(CC) $(CFLAGS) -o LAmerge LAmerge.c libalign.a -lm

LAshow: LAshow.c align.h DB.h QV.h
	$(CC) $(CFLAGS) -o LAshow LAshow.c libalign.a -lm

LAcat: LAcat.c align.h DB.h QV.h
	$(CC) $(CFLAGS) -o LAcat LAcat.c libalign.a -lm

LAsplit: LAsplit.c align.h DB.h QV.h
	$(CC) $(CFLAGS) -o LAsplit LAsplit.c libalign.a -lm

LAcheck: LAcheck.c align.h DB.h QV.h
	$(CC) $(CFLAGS) -o LAcheck LAcheck.c libalign.a -lm

LAupgrade.Dec.31.2014: LAupgrade.Dec.31.2014.c align.h DB.h QV.h
	$(CC) $(CFLAGS) -o LAupgrade.Dec.31.2014 LAupgrade.Dec.31.2014.c libalign.a -lm

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f LAupgrade.Dec.31.2014
	rm -f daligner.tar.gz
	rm -f libalign.a align.o DB.o QV.o
	rm -f libalign.so align.so DB.so QV.so

install: $(ALL)
	cp $(ALL) ~/bin

DESTDIR ?= ${HOME}/daligner

libinstall: libalign.a libalign.so
	mkdir -p ${DESTDIR}/include
	mkdir -p ${DESTDIR}/lib
	cp ${LIB_ALIGN_HEADERS} ${DESTDIR}/include/
	cp libalign.a libalign.so ${DESTDIR}/lib/
	chmod 755 ${DESTDIR}/lib/libalign.so

package:
	make clean
	tar -zcf daligner.tar.gz README *.h *.c Makefile

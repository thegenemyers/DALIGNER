CFLAGS = -O3 -Wall -Wextra -fno-strict-aliasing

ALL = daligner HPCdaligner HPCmapper LAsort LAmerge LAsplit LAcat LAshow LAcheck

all: $(ALL)

daligner: daligner.c filter.c filter.h align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o daligner daligner.c filter.c align.c DB.c QV.c -lpthread -lm

HPCdaligner: HPCdaligner.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o HPCdaligner HPCdaligner.c DB.c QV.c -lm

HPCmapper: HPCmapper.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o HPCmapper HPCmapper.c DB.c QV.c -lm

LAsort: LAsort.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAsort LAsort.c DB.c QV.c -lm

LAmerge: LAmerge.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAmerge LAmerge.c DB.c QV.c -lm

LAshow: LAshow.c align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAshow LAshow.c align.c DB.c QV.c -lm

LAcat: LAcat.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAcat LAcat.c DB.c QV.c -lm

LAsplit: LAsplit.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAsplit LAsplit.c DB.c QV.c -lm

LAcheck: LAcheck.c align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAcheck LAcheck.c align.c DB.c QV.c -lm

LAupgrade.Dec.31.2014: LAupgrade.Dec.31.2014.c align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAupgrade.Dec.31.2014 LAupgrade.Dec.31.2014.c align.c DB.c QV.c -lm

LIBRARY_HEADERS = align.h DB.h QV.h
LIBRARY_SOURCES = align.c DB.c QV.c
LIBRARY_OBJECTS = align.o DB.o QV.o
LIBRARY_SHAREDOBJECTS = align.so DB.so QV.so
SHARED_LIBRARY_SUFFIX=.so
SHARED_LIBRARY_EXTRA_LINKER_FLAGS=

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f LAupgrade.Dec.31.2014
	rm -f daligner.tar.gz
	rm -f ${LIBRARY_OBJECTS} ${LIBRARY_SHAREDOBJECTS} libalign.a libalign${SHARED_LIBRARY_SUFFIX}

install:
	cp $(ALL) ~/bin

package:
	make clean
	tar -zcf daligner.tar.gz README *.h *.c Makefile

libalign.a: ${LIBRARY_HEADERS} ${LIBRARY_SOURCES}
	gcc $(CFLAGS) -c -o align.o align.c
	gcc $(CFLAGS) -c -o DB.o DB.c
	gcc $(CFLAGS) -c -o QV.o QV.c
	ar -q libalign.a ${LIBRARY_OBJECTS}
	rm -f ${LIBRARY_OBJECTS}

libalign${SHARED_LIBRARY_SUFFIX}: ${LIBRARY_HEADERS} ${LIBRARY_SOURCES}
	gcc $(CFLAGS) -fpic -c -o align.so align.c
	gcc $(CFLAGS) -fpic -c -o DB.so DB.c
	gcc $(CFLAGS) -fpic -c -o QV.so QV.c
	gcc -shared ${SHARED_LIBRARY_EXTRA_LINKER_FLAGS} -o libalign${SHARED_LIBRARY_SUFFIX} ${LIBRARY_SHAREDOBJECTS}
	rm -f ${LIBRARY_SHAREDOBJECTS}

PREFIX=${HOME}
DESTDIR=

libinstall: libalign.a libalign${SHARED_LIBRARY_SUFFIX}
	mkdir -p ${DESTDIR}${PREFIX}/include
	cp ${LIBRARY_HEADERS} ${DESTDIR}${PREFIX}/include
	mkdir -p ${DESTDIR}${PREFIX}/lib
	cp libalign.a libalign${SHARED_LIBRARY_SUFFIX} ${DESTDIR}${PREFIX}/lib
	rm -f libalign.a libalign${SHARED_LIBRARY_SUFFIX}

CFLAGS = -O3 -Wall -Wextra -fno-strict-aliasing

ALL = daligner daligner_p HPCdaligner HPCmapper LAsort LAmerge LAsplit LAcat LAshow LAcheck LA4Falcon DB2Falcon DB.so

all: $(ALL)

daligner: daligner.c filter.c filter.h align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o daligner daligner.c filter.c align.c DB.c QV.c -lpthread -lm

daligner_p: daligner.c filter.c filter.h align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o daligner_p daligner.c filter.c align.c DB.c QV.c -lpthread -lm -DFALCON_DALIGNER_P

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

LA4Falcon: LA4Falcon.c align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LA4Falcon LA4Falcon.c align.c DB.c QV.c -lm

LAcat: LAcat.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAcat LAcat.c DB.c QV.c -lm

LAsplit: LAsplit.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAsplit LAsplit.c DB.c QV.c -lm

LAcheck: LAcheck.c align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAcheck LAcheck.c align.c DB.c QV.c -lm

DB2Falcon: DB2Falcon.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DB2Falcon DB2Falcon.c DB.c QV.c -lm

DB.so: DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -shared -fPIC -o DB.so DB.c QV.c -lm


clean:
	rm -f $(ALL)
	rm -f LAupgrade.Dec.31.2014
	rm -f daligner.tar.gz

LAupgrade.Dec.31.2014: LAupgrade.Dec.31.2014.c align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAupgrade.Dec.31.2014 LAupgrade.Dec.31.2014.c align.c DB.c QV.c -lm


install:
	cp $(ALL) ~/bin

package:
	make clean
	tar -zcf daligner.tar.gz README *.h *.c Makefile

CFLAGS = -O3 -Wall -Wextra -fno-strict-aliasing

all: daligner HPCdaligner \
     LAsort LAmerge LAsplit LAcat LAshow LAcheck

daligner: daligner.c filter.c filter.h align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o daligner daligner.c filter.c align.c DB.c QV.c -lpthread -lm

HPCdaligner: HPCdaligner.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o HPCdaligner HPCdaligner.c DB.c QV.c -lm

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

clean:
	rm -f daligner HPCdaligner
	rm -f LAsort LAmerge LAshow LAsplit LAcat LAcheck
	rm -f daligner.tar.gz

install:
	cp daligner HPCdaligner ~/bin
	cp LAsort LAmerge LAshow LAsplit LAcat LAcheck ~/bin

package:
	make clean
	tar -zcf daligner.tar.gz README *.h *.c Makefile

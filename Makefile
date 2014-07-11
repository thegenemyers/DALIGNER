CFLAGS = -O4 -Wall -Wextra

all: daligner HPCdaligner \
     LAsort LAmerge LAsplit LAcat LAshow LAcheck \

daligner: daligner.c filter.c filter.h align.c align.h DB.c DB.h
	gcc $(CFLAGS) -o daligner daligner.c filter.c align.c DB.c -lpthread -lm

HPCdaligner: HPCdaligner.c DB.c DB.h
	gcc $(CFLAGS) -o HPCdaligner HPCdaligner.c DB.c DB.h -lm

LAsort: LAsort.c align.h DB.c DB.h
	gcc $(CFLAGS) -o LAsort LAsort.c DB.c -lm

LAmerge: LAmerge.c align.h DB.c DB.h
	gcc $(CFLAGS) -o LAmerge LAmerge.c DB.c -lm

LAshow: LAshow.c align.c align.h DB.c DB.h
	gcc $(CFLAGS) -o LAshow LAshow.c align.c DB.c -lm

LAcat: LAcat.c align.h DB.c DB.h
	gcc $(CFLAGS) -o LAcat LAcat.c DB.c -lm

LAsplit: LAsplit.c align.h DB.c DB.h
	gcc $(CFLAGS) -o LAsplit LAsplit.c DB.c -lm

LAcheck: LAcheck.c align.c align.h DB.c DB.h
	gcc $(CFLAGS) -o LAcheck LAcheck.c align.c DB.c -lm

clean:
	rm -f daligner HPCdaligner
	rm -f LAsort LAmerge LAshow LAsplit LAcat LAcheck

install:
	cp daligner HPCdaligner ~/bin
	cp LAsort LAmerge LAshow LAsplit LAcat LAcheck ~/bin

package:
	make clean
	tar -zcf daligner.tar.gz README *.h *.c Makefile

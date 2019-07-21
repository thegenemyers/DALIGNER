DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = daligner HPC.daligner LAsort LAmerge LAsplit LAcat LAshow LAdump LAcheck LAa2b LAb2a dumpLA

all: $(ALL)

daligner: daligner.c filter.c filter.h radix.c radix.h align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o daligner -DWORD_SIZE=16  daligner.c filter.c radix.c align.c DB.c QV.c -lpthread -lm

HPC.daligner: HPC.daligner.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o HPC.daligner HPC.daligner.c DB.c QV.c -lm

LAsort: LAsort.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAsort LAsort.c DB.c QV.c -lm

LAmerge: LAmerge.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAmerge LAmerge.c DB.c QV.c -lm

LAshow: LAshow.c align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAshow LAshow.c align.c DB.c QV.c -lm

LAdump: LAdump.c align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAdump LAdump.c align.c DB.c QV.c -lm

LAcat: LAcat.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAcat LAcat.c DB.c QV.c -lm

LAsplit: LAsplit.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAsplit LAsplit.c DB.c QV.c -lm

LAcheck: LAcheck.c align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAcheck LAcheck.c align.c DB.c QV.c -lm

LAa2b: LAa2b.c align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAa2b LAa2b.c align.c DB.c QV.c -lm

LAb2a: LAb2a.c align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAb2a LAb2a.c align.c DB.c QV.c -lm

dumpLA: dumpLA.c align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o dumpLA dumpLA.c align.c DB.c QV.c -lm

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f daligner.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf daligner.tar.gz README.md Makefile *.h *.c

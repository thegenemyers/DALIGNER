DEST_DIR = ~/bin

CFLAGS += -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing
LDLIBS = -lm

ALL = daligner HPC.daligner LAsort LAmerge LAsplit LAcat LAshow LAdump LAcheck LAindex

all: $(ALL)

daligner: LDLIBS += -lpthread
daligner: daligner.c filter.c filter.h align.c align.h DB.c DB.h QV.c QV.h

HPC.daligner: HPC.daligner.c DB.c DB.h QV.c QV.h

LAsort: LAsort.c align.h DB.c DB.h QV.c QV.h

LAmerge: LAmerge.c align.h DB.c DB.h QV.c QV.h

LAshow: LAshow.c align.c align.h DB.c DB.h QV.c QV.h

LAdump: LAdump.c align.c align.h DB.c DB.h QV.c QV.h

LAcat: LAcat.c align.h DB.c DB.h QV.c QV.h

LAsplit: LAsplit.c align.h DB.c DB.h QV.c QV.h

LAcheck: LAcheck.c align.c align.h DB.c DB.h QV.c QV.h

LAupgrade.Dec.31.2014: LAupgrade.Dec.31.2014.c align.c align.h DB.c DB.h QV.c QV.h

LAindex: LAindex.c align.c align.h DB.c DB.h QV.c QV.h

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f LAupgrade.Dec.31.2014
	rm -f daligner.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf daligner.tar.gz README.md Makefile *.h *.c

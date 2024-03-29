DEST_DIR = ~/bin

# CFLAGS = -O0 -g -Wall -Wextra -Wno-unused-result -fno-strict-aliasing -fsanitize=address -fsanitize=undefined
# Above is for debug out of bound addresses, must compile with -lASAN -lUBSAN if gcc instead of clang

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = daligner HPC.daligner LAsort LAmerge LAsplit LAcat LAshow LA2ONE LAcheck ONE2LA

all: $(ALL)

daligner: daligner.c filter.c filter.h lsd.sort.c lsd.sort.h align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o daligner daligner.c filter.c lsd.sort.c align.c DB.c QV.c -lpthread -lm

HPC.daligner: HPC.daligner.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o HPC.daligner HPC.daligner.c DB.c QV.c -lm

LAsort: LAsort.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAsort LAsort.c DB.c QV.c -lm

LAmerge: LAmerge.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAmerge LAmerge.c DB.c QV.c -lm

LAshow: LAshow.c align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAshow LAshow.c align.c DB.c QV.c -lm

LA2ONE: LA2ONE.c align.c align.h DB.c DB.h QV.c QV.h ONElib.c ONElib.h
	gcc $(CFLAGS) -o LA2ONE LA2ONE.c align.c DB.c QV.c ONElib.c -lm

LAcat: LAcat.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAcat LAcat.c DB.c QV.c -lm

LAsplit: LAsplit.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAsplit LAsplit.c DB.c QV.c -lm

LAcheck: LAcheck.c align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAcheck LAcheck.c align.c DB.c QV.c -lm

ONE2LA: ONE2LA.c align.c align.h DB.c DB.h QV.c QV.h ONElib.c ONElib.h
	gcc $(CFLAGS) -o ONE2LA ONE2LA.c align.c DB.c QV.c ONElib.c -lm

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f daligner.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf daligner.tar.gz README.md Makefile *.h *.c


MANPAGES = $(basename $(wildcard  *.1.md))

all: $(MANPAGES)

%.1: %.1.md
	pandoc -s -f markdown -t man $< > $@

clean:
	$(RM) $(MANPAGES)

.PHONY: clean

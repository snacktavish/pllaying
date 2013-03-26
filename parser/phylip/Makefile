CC           = gcc
CFLAGS       = -g -Wall 

all: parser

parser: phylip.c ssort.c tester.c ../../lexer.c
	@echo "===> COMPILE $@"
	$(CC) $(CFLAGS) $+ -o $@

clean:
	-rm parser

.PHONY: all clean

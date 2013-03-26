CC           = gcc
CFLAGS       = -g -Wall

all: parser

parser: newick.c tester.c ../../lexer.c ../../mem_alloc.c ../../stack.c
	@echo "===> COMPILE $@"
	$(CC) $(CFLAGS) $+ -o $@ -lm

clean:
	-rm parser

.PHONY: all clean

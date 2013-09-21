CC           = gcc
CFLAGS       = -g -Wall -D__DEBUGGING_MODE

all: parser

parser: newick.c tester.c ../common.c ../../lexer.c ../../stack.c
	@echo "===> COMPILE $@"
	$(CC) $(CFLAGS) $+ -o $@ -lm

clean:
	-rm parser

.PHONY: all clean

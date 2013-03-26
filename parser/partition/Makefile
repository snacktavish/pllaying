CC           = gcc
CFLAGS       = -g -Wall

all: parser

parser: part.c tester.c ../../lexer.c ../../mem_alloc.c ../../queue.c ../../hash.c
	@echo "===> COMPILE $@"
	$(CC) $(CFLAGS) $+ -o $@ -lm

clean:
	-rm parser

.PHONY: all clean

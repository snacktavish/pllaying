CC           = gcc
CFLAGS       = -g -Wall 

all: parser

parser: lexer.c phylip.c ssort.c main.c
	@echo "===> COMPILE $@"
	$(CC) $(CFLAGS) $+ -o $@

clean:
	-rm parser

.PHONY: all clean

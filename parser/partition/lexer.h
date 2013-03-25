#ifndef LEXER_H
#define LEXER_H

#define  SIZE_ASCII              128
#define  EOS                     0x00000200

#define  SYMBOL_CR               1 << 0
#define  SYMBOL_LF               1 << 1
#define  SYMBOL_LFCR             1 << 2
#define  SYMBOL_DIGIT            1 << 3
#define  SYMBOL_CHAR             1 << 4
#define  SYMBOL_SPACE            1 << 5
#define  SYMBOL_TAB              1 << 6
#define  SYMBOL_EOF              1 << 7
#define  SYMBOL_UNKNOWN          1 << 8
#define  SYMBOL_DOT              1 << 9
#define  SYMBOL_COLON            1 << 10
#define  SYMBOL_OPAREN           1 << 11
#define  SYMBOL_CPAREN           1 << 12
#define  SYMBOL_COMMA            1 << 13
#define  SYMBOL_SEMICOLON        1 << 14
#define  SYMBOL_EQUAL            1 << 15
#define  SYMBOL_DASH             1 << 16
#define  SYMBOL_SLASH            1 << 17

#define  LEX_NUMBER              1 << 0
#define  LEX_STRING              1 << 1
#define  LEX_EOF                 1 << 2
#define  LEX_WHITESPACE          1 << 3
#define  LEX_NEWLINE             1 << 4
#define  LEX_UNKNOWN             1 << 5
#define  LEX_COLON               1 << 6
#define  LEX_OPAREN              1 << 7
#define  LEX_CPAREN              1 << 8
#define  LEX_FLOAT               1 << 9
#define  LEX_COMMA               1 << 10
#define  LEX_SEMICOLON           1 << 11
#define  LEX_EQUAL               1 << 12
#define  LEX_DASH                1 << 13
#define  LEX_SLASH               1 << 14

struct ltoken_t
 {
   int          class;
   const char * lexeme;
   int          len;
 };

int get_next_byte (void);
int get_next_symbol (void);
struct ltoken_t get_token (int * input);
void init_lexan (const char * text, int n);

#endif

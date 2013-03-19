#ifndef LEXER_H
#define LEXER_H

#define  SIZE_ASCII              128
#define  EOS                     0x00000200

#define  SYMBOL_CR               0x00000001
#define  SYMBOL_LF               0x00000002
#define  SYMBOL_LFCR             0x00000004
#define  SYMBOL_DIGIT            0x00000008
#define  SYMBOL_CHAR             0x00000010
#define  SYMBOL_SPACE            0x00000020
#define  SYMBOL_TAB              0x00000040
#define  SYMBOL_EOF              0x00000080
#define  SYMBOL_UNKNOWN          0x00000100
#define  SYMBOL_DOT              0x00000200
#define  SYMBOL_COLON            0x00000400
#define  SYMBOL_OPAREN           0x00000800
#define  SYMBOL_CPAREN           0x00001000
#define  SYMBOL_COMMA            0x00002000
#define  SYMBOL_SEMICOLON        0x00004000

#define  LEX_NUMBER              0x00000001
#define  LEX_STRING              0x00000002
#define  LEX_EOF                 0x00000004
#define  LEX_WHITESPACE          0x00000008
#define  LEX_NEWLINE             0x00000010
#define  LEX_UNKNOWN             0x00000020
#define  LEX_COLON               0x00000040
#define  LEX_OPAREN              0x00000080
#define  LEX_CPAREN              0x00000100
#define  LEX_FLOAT               0x00000200
#define  LEX_COMMA               0x00000400
#define  LEX_SEMICOLON           0x00000800

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

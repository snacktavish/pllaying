#include <stdio.h>
#include "lexer.h"

static const char * rawtext;
static int          rawtext_size;
static int          pos = 0;

int lex_table[SIZE_ASCII] = {
/*      */ SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/*      */ SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/*      */ SYMBOL_UNKNOWN,     SYMBOL_TAB,      SYMBOL_CR, SYMBOL_UNKNOWN,
/*      */ SYMBOL_UNKNOWN,      SYMBOL_LF, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/*      */ SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/*      */ SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/*      */ SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/*      */ SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/*  !"# */   SYMBOL_SPACE, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/* $%&' */ SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/* ()*+ */  SYMBOL_OPAREN,  SYMBOL_CPAREN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/* ,-./ */   SYMBOL_COMMA,    SYMBOL_DASH,     SYMBOL_DOT,   SYMBOL_SLASH,
/* 0123 */   SYMBOL_DIGIT,   SYMBOL_DIGIT,   SYMBOL_DIGIT,   SYMBOL_DIGIT,
/* 4567 */   SYMBOL_DIGIT,   SYMBOL_DIGIT,   SYMBOL_DIGIT,   SYMBOL_DIGIT,
/* 89:; */   SYMBOL_DIGIT,   SYMBOL_DIGIT,   SYMBOL_COLON, SYMBOL_SEMICOLON,
/* <=>? */ SYMBOL_UNKNOWN,   SYMBOL_EQUAL, SYMBOL_UNKNOWN,    SYMBOL_CHAR,
/* @ABC */ SYMBOL_UNKNOWN,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* DEFG */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* HIJK */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* LMNO */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* PQRS */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* TUVW */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* XYZ[ */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR, SYMBOL_UNKNOWN,
/* \]^_ */ SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,    SYMBOL_CHAR,
/* `abc */ SYMBOL_UNKNOWN,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* defg */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* hijk */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* lmno */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* pqrs */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* tuvw */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* xyz{ */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR, SYMBOL_UNKNOWN,
/* |}~  */    SYMBOL_CHAR, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN
 };

extern int debug ;

int 
get_next_byte (void)
{
  if (pos == rawtext_size) 
   {
     ++pos;
     return (EOS);
   }

  return (rawtext[pos++]);
}

int
get_next_symbol ()
{
  int ch, sym;

  ch = get_next_byte ();

  if (ch == EOS) return (SYMBOL_EOF);
  if (ch >= SIZE_ASCII) return SYMBOL_UNKNOWN;

  sym = lex_table[ch];

  if (sym == SYMBOL_LF)
   {
     if (get_next_byte() == '\n')
      {
        sym = SYMBOL_LFCR;
      }
     else
      {
       --pos;
      }
   }

  return sym;
}

struct ltoken_t
get_token (int * input)
{
  struct ltoken_t   token;
  int               start_pos;
  int               floating = 0;

  token.lexeme = rawtext + pos - 1;
  start_pos    = pos;

  switch (*input)
   {
     case SYMBOL_SLASH:
       token.class = LEX_SLASH;
       *input = get_next_symbol();
       break;

     case SYMBOL_DASH:
       token.class = LEX_DASH;
       *input = get_next_symbol();
       break;

     case SYMBOL_EQUAL:
       token.class = LEX_EQUAL;
       *input = get_next_symbol();
       break;

     case SYMBOL_SEMICOLON:
       token.class = LEX_SEMICOLON;
       *input = get_next_symbol();
       break;

     case SYMBOL_COMMA:
       token.class = LEX_COMMA;
       *input = get_next_symbol();
       break;

     case SYMBOL_COLON:
       token.class = LEX_COLON;
       *input = get_next_symbol();
       break;

     case SYMBOL_OPAREN:
       token.class = LEX_OPAREN;
       *input = get_next_symbol();
       break;

     case SYMBOL_CPAREN:
       token.class = LEX_CPAREN;
       *input = get_next_symbol();
       break;

     case SYMBOL_SPACE:
     case SYMBOL_TAB:
       do
        {
          *input = get_next_symbol();
        } while (*input == SYMBOL_SPACE || *input == SYMBOL_TAB);
       token.len   = pos - start_pos;
       token.class = LEX_WHITESPACE; 
       if (*input == SYMBOL_LFCR) --token.len;
       break;
       
     case SYMBOL_DIGIT:
       do
        {
          *input = get_next_symbol();   
        } while (*input == SYMBOL_DIGIT);

       if (*input == SYMBOL_DOT)
        {
          floating = 1;
          do
           {
             *input = get_next_symbol ();
           } while (*input == SYMBOL_DIGIT);
        }

       if (*input != SYMBOL_CHAR)
        {
          token.len   = pos - start_pos;
          if (!floating)
            token.class = LEX_NUMBER;
          else
            token.class = LEX_FLOAT;
        }
       else
        {
          do {
            *input = get_next_symbol();
          } while (*input == SYMBOL_CHAR || *input == SYMBOL_DIGIT || *input == SYMBOL_DOT);
          token.len   = pos - start_pos;
          token.class = LEX_STRING;
        }

       if (*input == SYMBOL_LFCR) --token.len;
       break;

     case SYMBOL_CHAR:
       do
        {
          *input = get_next_symbol();
        } while (*input == SYMBOL_CHAR || *input == SYMBOL_DIGIT || *input == SYMBOL_DASH);
       token.len   = pos - start_pos;
       token.class = LEX_STRING;
       if (*input == SYMBOL_LFCR) --token.len;
       break;
       
     case SYMBOL_EOF:
       token.class = LEX_EOF;
       break;

     case SYMBOL_CR:
     case SYMBOL_LF:
     case SYMBOL_LFCR:
       do
        {
          *input = get_next_symbol();
        } while (*input == SYMBOL_CR || *input == SYMBOL_LFCR || *input == SYMBOL_LF);
       token.class = LEX_NEWLINE;
       break;
     case SYMBOL_UNKNOWN:
       token.class = LEX_UNKNOWN;
       break;
   }

  return (token);
}

void
init_lexan (const char * text, int n)
{
  rawtext      = text;
  rawtext_size = n;
  pos          = 0;
}
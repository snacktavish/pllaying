#include <stdio.h>
#include "../mem_alloc.h"
#include "common.h"

/** @brief Read the contents of a file
    
    Reads the ile \a filename and return its content. In addition
    the size of the file is stored in the input variable \a filesize.
    The content of the variable \a filesize can be anything and will
    be overwritten.

    @param filename
      Name of the input file

    @param filesize
      Input parameter where the size of the file (in bytes) will be stored

    @return
      Contents of the file
*/
char * 
pllReadFile (const char * filename, long * filesize)
{
  FILE * fp;
  char * rawdata;

  fp = fopen (filename, "r");
  if (!fp) return (NULL);

  /* obtain file size */
  if (fseek (fp, 0, SEEK_END) == -1)
   {
     fclose (fp);
     return (NULL);
   }

  *filesize = ftell (fp);

  if (*filesize == -1) 
   {
     fclose (fp);
     return (NULL);
   }
  rewind (fp);

  /* allocate buffer and read file contents */
  rawdata = (char *) rax_malloc (((*filesize) + 1) * sizeof (char));
  if (rawdata) 
   {
     if (fread (rawdata, sizeof (char), *filesize, fp) != (size_t) *filesize) 
      {
        rax_free (rawdata);
        rawdata = NULL;
      }
     else
      {
        rawdata[*filesize] = 0;
      }
   }

  fclose (fp);

  return (rawdata);
}


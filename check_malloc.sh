#!/bin/bash

# check if the executable wants to dynamically import any of the malloc/free
# related symbols.
 
nm -D $1 | egrep -e "U malloc$" -e "U memalign$" -e "U posix_memalign$" -e "U free$" -e "U realloc$" -e "U calloc$"

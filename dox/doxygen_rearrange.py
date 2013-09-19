#!/usr/bin/python

import re
import sys


g_memberdecls_table_string = ('table', 'class', 'memberdecls')
g_bStatic = 0
g_bNormal = 1
g_static_funtable_start = '<table class="memberdecls"><tr class="heading"><td colspan="2"><h2 class="groupheader"><a name="func-members"></a>Static functions</h2></td></tr>\n'
g_static_funtable_end   = '</table>\n'

def locate_memberdecls_index (content):
  dfa = re.compile (".*".join (g_memberdecls_table_string))
  for i in range(len(content)):
    if dfa.search (content[i]) is not None: break

  return i
   

def readFile (fname):
  with open(fname) as f:
    content = f.readlines()

  return content

def split_html_content (content):
  i = locate_memberdecls_index (content)

  # nothing found or found in the last line, quit
  if i == len (content) - 1: return (content, [], [])

  dfa = re.compile ("^</table>")
  # get the ending </table> tag
  for j in range (i + 3, len(content)):
    if dfa.search (content[j]) is not None: break
  
  header    = content[:i + 3]
  functions = content[i + 3:j]
  footer    = content[j:]

  return (header, functions, footer)

def split_static_functions (functions):
  
  dfa_fun  = re.compile ("^<tr[^>]*class[^>]*memitem[^>]*><td[^>]*>static")
  dfa_desc = re.compile ("^<tr[^>]*class[^>]*memdesc[^>]*>")
  static = []
  normal = []

  where = 0

  i = 0
  while (i < len(functions)):
    while len (functions[i].strip()) == 0: i = i + 1
    if dfa_fun.search (functions[i]) is not None:
      where = g_bStatic
      static.append (functions[i])
    else:
      normal.append (functions[i])
      where = g_bNormal

    i = i + 1

    while len (functions[i].strip()) == 0: i = i + 1
    # check if there is a @brief description and append it to the appropriate list
    if dfa_desc.search (functions[i]) is not None:
      if where == g_bStatic:
        static.append (functions[i])
      else:
        normal.append (functions[i])
      i = i + 1

    while len (functions[i].strip()) == 0: i = i + 1
    # we assume there is a horizontal bar
    if where == g_bStatic:
      static.append (functions[i])
    else:
      normal.append (functions[i])

    i = i + 1
  
  return static, normal

def overwrite_file (fname, header, static, normal, footer):
  
  f = open (fname, 'w')

  normalExist = True

  for line in header:
    f.write(line)


  if len(normal) != 0:
    for line in normal:
      f.write(line)
  else:
    normalExist = False

  if normalExist == True and len(static) != 0:
    f.write(g_static_funtable_start)
    for line in static:
      f.write(line)
    f.write(g_static_funtable_end)
  else:
    if normalExist == False:
      for line in static:
        f.write(line)

  for line in footer:
    f.write(line)
  
  f.close()


def process_static_functions (fname):
  """Separate static from normal functions in file fname and overwrite the file"""
  content = readFile (fname)
  header, functions, footer = split_html_content (content)

  # no functions found? quit
  if not functions: return
  
  static, normal = split_static_functions (functions)

  overwrite_file (fname, header, static, normal, footer)


if __name__ == "__main__":
  process_static_functions (sys.argv[1])

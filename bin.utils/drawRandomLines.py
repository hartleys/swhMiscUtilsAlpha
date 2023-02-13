#!/usr/bin/env python
from __future__ import print_function

import sys,gzip,errno,time,random,collections,re,random

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    return False

def intOrZero(s):
    try:
        return int(s)
    except ValueError:
        pass
    return 0

def floatOrZero(s):
    try:
        return float(s)
    except ValueError:
        pass
    return 0

def openSmart(s,t):
   if(s[-3:] == ".gz"):
      return gzip.open(s,t);
   else:
      return open(s,t);

def indexOf(lst,x):
    try:
       y = lst.index(x);
       return y;
    except ValueError:
       return -1;

def dequote(s):
    #If a string has single or double quotes around it, remove them.
    #Make sure the pair of quotes match.
    #If a matching pair of quotes is not found, return the string unchanged.
    if (s[0] == s[-1]) and s.startswith(("'", '"')):
        return s[1:-1]
    return s


eprint("STARTING");

args = sys.argv;

#+SHORTDESC="Takes an input stream, spits out a subset of lines, in the input order."
#+ HELPTEXT="..."$'\n'
#++HELPTEXT=" "$'\n'
#+ SYNTAXTEXT="drawRandomLines.py <fraction>"$'\n'
#+ PARAMS="--help: displays syntax and help info."$'\n'
#++PARAMS="--man: synonym for --help."$'\n'
#+ VERSION="0.0.5"

if len(args) == 0 or "--help" in args or "--man" in args:
   import os;
   scriptAbsPath=os.path.abspath(__file__);
   shareMiscHome=os.path.dirname(os.path.abspath( os.path.dirname( __file__ )));
   os.system(shareMiscHome+"/internal.manUtil/manForScript "+scriptAbsPath);
   exit(1);

keepLines = 0;
if "--keepLines" in args:
   idx = args.index("--keepLines");
   args.pop(idx);
   keepLines = int(args.pop(idx));

eprint("INDATA="+args[-1]);


inf = sys.stdin;

pct = float(args[-1]);
out = sys.stdout;

lnct = 0;
keepct = 0;
for line in inf:
    lnct = lnct + 1;
    if lnct % 1000 == 0:
      eprint(".", end='');
      if lnct % 5000 == 0:
         eprint(" ",end='');
      if lnct % 10000 == 0:
         cells = line.split("\t")
         eprint("["+str(lnct)+" lines complete] [kept:",keepct,"/",lnct,"] ["+time.strftime("%c")+"]")
    if random.random() < pct or lnct <= keepLines:
       out.write(line);
       keepct = keepct + 1;


eprint("["+str(lnct)+" lines complete] [kept:",keepct,"/",lnct,"] ["+time.strftime("%c")+"]")


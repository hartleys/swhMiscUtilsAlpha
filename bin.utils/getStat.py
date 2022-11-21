#!/usr/bin/env python
from __future__ import print_function

import sys,gzip,errno,time

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    return False

args = sys.argv;

#+SHORTDESC="Gets basic stats from a character stream"
#+ HELPTEXT="Get statistics from stream (min, max, sum, mean)"$'\n'
#++HELPTEXT=" "$'\n'
#+ SYNTAXTEXT="getStat.py [options]"$'\n'
#+ PARAMS="--help: displays syntax and help info."$'\n'
#++PARAMS="--man: synonym for --help."$'\n'
#+ VERSION="0.0.5"

getMax = False
getMin = False
getMean = False
getSum = False
if "--getMax" in args:
   getMax = True;
if "--getMin" in args:
   getMin = True;
if "--getMean" in args:
   getMean = True;
if "--getSum" in args:
   getSum = True;

asFloat = False;
if "--float" in args:
   asFloat = True;

inf = sys.stdin;
out = sys.stdout;

lnct = 0;
currMin = -1
currMax = -1
currSum = 0


try:
   for line in inf:
      currString = line[:-1];
      if isNumber(currString):
         curr = float(currString)
         if lnct == 0:
            currMin = curr
            currMax = curr
         else:
            currMin = min(curr,currMin);
            currMax = max(curr,currMax);
         currSum = currSum + curr;
         lnct = lnct + 1;
      else:
         eprint("Warning: non-numeric found: \""+currString+"\"")
except IOError as e:
   if e.errno == errno.EPIPE:
      exit(0);

outList = []

if asFloat:
  if getMax:
     outList = outList + [currMax]
  if getMin:
     outList = outList + [currMin]
  if getSum:
     outList = outList + [currSum]
  if getMean:
     outList = outList + [float(currSum) / float(lnct)]
else:
  if getMax:
     outList = outList + [int(currMax)]
  if getMin:
     outList = outList + [int(currMin)]
  if getSum:
     outList = outList + [int(currSum)]
  if getMean:
     outList = outList + [int(float(currSum) / float(lnct))]

out.write("\t".join([str(s) for s in outList]));




















#!/bin/bash

#estat | tail -n+2 | awk '{ if ( $5 == "qw" ) print $1 }' | tr '\n' '|'

#TODO-WRITE-HELPDOC
#+SHORTDESC="TODO-WRITE-HELPDOC"
#+ HELPTEXT="..."$'\n'
#++HELPTEXT="More Help Text"$'\n'
#+ SYNTAXTEXT="... [options]"$'\n'
#+ PARAMS="--help: displays syntax and help info."$'\n'
#++PARAMS="--man: synonym for --help."$'\n'
#+ VERSION="0.0.6"
CURRENT_SCRIPT_DIR="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd -P )"
if [ "$1" == "--help" -o "$1" == "--man" ]; then
  ${CURRENT_SCRIPT_DIR}/../internal.manUtil/manForScript ${BASH_SOURCE[0]}
  exit 1
elif [ "$1" == "--manMarkdown" ]; then
  ${CURRENT_SCRIPT_DIR}/../internal.manUtil/manForScript --markdown ${BASH_SOURCE[0]}
  exit 1
fi


if [ "$SWHENV_QUEUEENGINE_QUEUED" == "" ]; then
   SWHENV_QUEUEENGINE_QUEUED="qw";
fi


##PYTHON:
if "--help" in args or "--man" in args:
   import inspect,commands,os
   SCRIPT_FILE_PATH=inspect.stack()[0][1];
   print(SCRIPT_FILE_PATH);
   SCRIPT_DIR_PATH=commands.getstatusoutput("dirname "+SCRIPT_FILE_PATH)[1]
   SCRIPT_FILENAME=commands.getstatusoutput("basename "+SCRIPT_FILE_PATH)[1]
   #print("SCRIPT_DIR_PATH="+SCRIPT_DIR_PATH);
   #print("SCRIPT_FILENAME="+SCRIPT_FILENAME);
   os.system(SCRIPT_DIR_PATH+"/../internal.manUtil/manForScript"+" "+SCRIPT_FILE_PATH)
   exit(0);

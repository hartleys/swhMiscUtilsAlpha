#!/bin/bash

CURRENT_SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

while read line
do
  if [ "$(head -n1 $line | cut -c1-2)" == '#!' ]; then
    SHORTDESC="";
    echo "line=$line" >2
    source <( cat $line | grep "^#[+]" | sed -e 's/^#[+]//' | sed -r 's/[+]([A-Za-z_]+)="/\1="${\1}/'   );
    if [ "$SHORTDESC" == "" ]; then
      echo "$(basename $line)"$'\t'"(no description)"
    else
      echo "$(basename $line)"$'\t'"$SHORTDESC"
    fi
  fi
done < <(find $CURRENT_SCRIPT_DIR/../* -maxdepth 0 -perm -u+x -type f) | column -s$'\t' -t > $CURRENT_SCRIPT_DIR/COMMANDLIST.txt

cat $CURRENT_SCRIPT_DIR/README_START > $CURRENT_SCRIPT_DIR/README
echo "" >> $CURRENT_SCRIPT_DIR/README
echo "## Available Commands:" >> $CURRENT_SCRIPT_DIR/README
cat $CURRENT_SCRIPT_DIR/COMMANDLIST.txt >> $CURRENT_SCRIPT_DIR/README
echo "" >> $CURRENT_SCRIPT_DIR/README
cat $CURRENT_SCRIPT_DIR/README_END >> $CURRENT_SCRIPT_DIR/README



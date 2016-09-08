#!/bin/bash

set -e

mode=""
file=""
tmpdir=tmptmp
fixed=false


function clean_up {
   if [[ $? == 0 ]]; then
      echo "[happy] Script succeeded"
   else
      echo "[grumpy] Script failed"
   fi
   if [[ $tmpdir =~ tmptmp ]]; then
      rm -rf $tmpdir
      echo "I've removed the working directory $tmpdir"
   else
      echo "I've left the working directory $tmpdir in place"
   fi
}


   ##  register the function clean_up to be run when script exits,
   ##  after all the other code in this script.
   ##
trap clean_up SIGTERM  EXIT


   ## This is some boilerplate code to parse arguments. Options
   ## are all single-character.
   ## Options followed by a colon expect an argument;
   ## I've grouped all the (two) unary options (those that do not expect an argument)
   ## at the end.
   ## The leading colon helps with error processing using :) and ?) below.
   ## Just use this as a boilerplate example.
   ##
while getopts :f:m:t:Fh opt
do
    case "$opt" in
    m)
      mode=$OPTARG
      ;;
    f)
      file=$OPTARG
      ;;
    t)
      tmpdir=$OPTARG
      ;;
    F)
      fixed=true
      ;;
    h)
      cat <<EOU
-f file
-m mode
EOU
      exit
      ;;
    :) echo "Flag $opt needs argument"
        exit 1;;
    ?) echo "Flag $opt unknown"
        exit 1;;
   esac
done


if [[  $file == "" ]]; then
   echo "Need file! (see -h)"
   false
fi

if $fixed; then
   echo "Running in fixed mode"
fi

OPTIND=$(($OPTIND-1))
shift $OPTIND


for arg in $@; do
   echo $arg
done


mkdir -p $tmpdir

cp "$file" $tmpdir

nl=$(wc -l < $tmpdir/$file)

echo "file $file has $nl lines"




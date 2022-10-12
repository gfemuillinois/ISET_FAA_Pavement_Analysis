#!/bin/bash

# STEP 0) Set place to look for dynamic libraries
export LD_LIBRARY_PATH=`pwd`/lib # this has to be first otherwise the user computer might search wrong for ldd

# STEP 1) Check and change .so libnames to match the ones the executable expects

#function to looks for a string on a table of strings
function lookstring {
  isFound=false
  for (( is=2; is<=$#; is++))
  do
    num=$is
    istring=${!num}
    if [ "$istring" = "$1" ]; then
      isFound=true
      break
    fi
  done
}

for (( it=1; it<=2; it++)) # two times to guarantee
do
  # Getting all files in folder lib
  FILES=`pwd`/lib/*

  # Executing command ldd on executable to know the names expected from tcliset
  MYLDD=$(ldd `pwd`/executable/tcliset)
  unset LDDONLYLIB

  # Getting only the lib names from "ldd tcliset" command
  for ldd in $MYLDD
  do
    if [[ $ldd == *"lib"* && $ldd != *"linux"* && $ldd != *"usr"* && $ldd != *"ISET"* ]]
    then
      LDDONLYLIB="$LDDONLYLIB $ldd"
    fi
  done

  # Going through all the libs and changing the name in case they are different
  # from what is expected by tcliset
  for f in $FILES
  do
    remainingname="${f%.*}"
    myext="${f##*.}"
    if [ "$myext" = "so" ]; then
      continue
    else
      remainingname2="${remainingname%.*}"
      myext2="${remainingname##*.}"
      filename="${remainingname##*/}"
      lookstring $filename $LDDONLYLIB
      if [ $isFound = true ]; then
        #echo "found $remainingname"
        mv $f $remainingname
        continue
      fi
      if [ "$myext2" = "so" ]; then
        continue
      else
        remainingname3="${remainingname2%.*}"
        myext3="${remainingname2##*.}"      
        filename="${remainingname2##*/}"
        lookstring $filename $LDDONLYLIB
        if [ $isFound = true ]; then
          #echo "found2 $remainingname2"
          mv $f $remainingname2
          continue
        fi
        if [ "$myext3" = "so" ]; then # evertything after here is useless for now. It was used for debugging (consider erasing)
          #echo "found so..."
          continue
        else
          remainingname4="${remainingname3%.*}"
          myext4="${remainingname3##*.}" 
          if [ "$myext4" = "so" ]; then
            #echo "found so 2..."
            continue
          fi  
        fi
      fi
    fi
  done #end for FILES
done #end for it
# STEP 2) Printing ISET information

# Clear screen
clear

# Set prompt
export PS1="\[\033[01;32m\]\u@\h:\[\033[01;34m\]\W\$ \[\033[00m\]"

# Set terminal title
echo -n -e "\033]0;ISET - The Ultimate GFEM Library \007"

cat << EOF
        _____ 
       |_   _| _____       _______
         | |  |  ___| ___ |__   __|
         | |  | |___ / _ \   | |
        _| |_ |___  |  __/   | |
        \___/ |_____|\___|   |_| 
    ==================================

This is a shell with PATHs and EXTERNAL_LIBs setup to work with ISET.
To run iset type anywhere where there is a tcl 

    % iset NameOfTcl.tcl

For a first test, go to folder examples/1_FirstTest and type

    % iset 2d_tria_heat_dirich_pt_analytic_p1.tcl    

For more information please send an email to one of following:

Nathan Shauer               - nathan.sh@gmail.com
Prof. Carlos Armando Duarte - caduarte2@gmail.com

EOF

# STEP 3) Creating alias to execute iset anywhere
alias iset=`pwd`/executable/tcliset


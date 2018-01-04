#!/bin/bash

  # bail is a simple error exit routine
bail(){
    echo
    echo "================"
    echo "      ERROR:    "
    echo "$*"   
    echo "    END ERROR   "
    echo "================"
    echo
    
    exit 1
  }

acc_cp() {
    # Attempt to access a file and if that fails attempt to copy it
    [ -z "$1" ] && bail "--> acc_cp <-- requires a local file name as the first arg"
    local_file=$1
    [ -z "$2" ] && bail "--> acc_cp $1 <-- requires a remote file name as the second arg"
    remote_file=$2
    shift 2
    # Alternate path to a directory where .queue/.crawork will be found
    JHOME=''

    access $local_file $remote_file $* ||
      cp -f $remote_file $local_file ||
      cp -f ${INPUT_PATH}$remote_file $local_file ||
      bail "--> acc_cp $local_file $remote_file $* <-- Unable to locate $remote_file"

    sha1sum $local_file >> xnemo_access_log
  }

  # mod_nl is used to change namelist values for the current job
mod_nl(){
    # Modify a text file containing namelists
    # Each variable definition in this input namelist file must be on a line by itself
    # and be the of the form var =.* which will be replaced with var = value
    #
    # Usage: mod_nl namelist_file_in var1 [var2 ...]
    #   namelist_file_in  ...is an existing text file containing the namelists
    #   var1 var2 ...     ...names of variables to be changed
    #                        these variables must be defined in the current environment
    #
    [ -z "$1" ] && bail "mod_nl requires a namelist file as the first arg"
    [ -z "$2" ] && bail "mod_nl requires at least 1 variable name on the command line"
    nlfile=$1
    [ ! -s $nlfile ] && bail "mod_nl input namelist file --> $nlfile <-- is missing or empty"
    shift

    # Create a backup copy of the input namelist file, overwriting any existing backup
    cp -f $nlfile ${nlfile}.bak

    # Write a sed program to make the requested substitutions
    sprog=mod_nl_cmd
    rm -f $sprog
    touch $sprog
    for var in $*; do
      eval val=\$$var
      [ -z "$val" ] && bail "mod_nl $var is not defined or null"
      if [ x`echo $var|sed -n '/^ *c/p'` != x ]; then
        # This is a character variable that needs to be quoted
        # This assumes the nemo naming convention (char variables start with "c")
        vsub='s/^ *'$var' *=.*/'"  ${var} = \"$val\""'/'
      else
        vsub='s/^ *'$var' *=.*/'"  ${var} = $val"'/'
      fi
      # Add this command to the file
      # Note, the quotes are required to preserve white space
      echo "$vsub" >> $sprog

      # Echo the changes to stdout (useful for debugging runs)
      echo "${nlfile}: RESET $var = $val"
    done

    # Make the substitutions in situ, overwriting the input file
    sed -f $sprog $nlfile > mod_nl_tmp
    mv mod_nl_tmp $nlfile
    rm $sprog
  }  

  ToF(){
    #   usage: ToF var_name
    # purpose: Possibly reset the value of var_name to "0" (false) or "1" (true)
    #          If var_name is null or has a value of "off" or "no" then reset to "0"
    #          If var_name has a value of "on" or "yes" then reset to "1"
    #          Otherwise return with var_name unchanged
    [ -z "$1" ] && bail "ToF requires a variable name as an argument"
    eval ToF_var\=\$$1
    XXX=`echo $ToF_var|sed 's/ //g'`
    eval ToF_var\=$XXX
    if [ -n "$ToF_var" ]; then
      if   [ "$ToF_var" = 'on'  ]; then eval ToF_var\=1
      elif [ "$ToF_var" = 'off' ]; then eval ToF_var\=0
      elif [ "$ToF_var" = 'yes' ]; then eval ToF_var\=1
      elif [ "$ToF_var" = 'no'  ]; then eval ToF_var\=0
      else
        eval ToF_var\=\$$1
      fi
    else
      eval ToF_var\=0
    fi
    eval $1=$ToF_var
  }

release(){
    rm -f $1
  }

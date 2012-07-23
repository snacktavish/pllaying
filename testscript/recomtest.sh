#! /bin/sh
# This is the working directory you want to run the test from
WORKINGDIR=`pwd`
DATADIR=`pwd`/../testdata
# Set the number of threads you want to use
NUMPTHREADS=4
# NO NEED TO EDIT BEYOND THIS POINT
SSE3_GCC="testscript.SSE3.gcc"
SSE3_PTHREADS_GCC="testscript.SSE3.PTHREADS.gcc"

# KEEP A LOG OF STDOUT OF RAxML
LOGFILE=RAxML_log.txt
ERRLOGFILE=RAxML_ERR_log.txt

run_SSE3_GCC()
{
  sh ${WORKINGDIR}/${SSE3_GCC} ${OPTIONS}
  RETVAL=$?
  if [ $RETVAL -ne 0 ] ; then
    echo " Error: testscript.SSE3.gcc failed. Check the messages above ...";
    echo "\n"
    echo " Supertestscript did not finish successfully."
    echo "\n"
    exit 1
  fi
  echo "\n"
}

run_SSE3_PTHREADS_GCC()
{
  sh ${WORKINGDIR}/${SSE3_PTHREADS_GCC} ${OPTIONS} -T ${NUMPTHREADS}
  RETVAL=$?
  if [ $RETVAL -ne 0 ] ; then
    echo " Error: testscript.SSE3.PTHREADS.gcc failed. Check the messages above ...";
    echo "\n"
    echo " Supertestscript did not finish successfully."
    echo "\n"
    exit 1
  fi
  echo "\n"
}


if [ ! -f ${WORKINGDIR}/${SSE3_GCC} ] ; then
  echo " Error: file ${SSE3_GCC} is missing from ${WORKINGDIR}, exiting ...";
  exit 1
fi
if [ ! -f ${WORKINGDIR}/${SSE3_PTHREADS_GCC} ] ; then
  echo " Error: file ${SSE3_PTHREADS_GCC} is missing from ${WORKINGDIR}, exiting ...";
  exit 1
fi

if [ $# -eq 0 ] ; then
  echo "\n"
  echo " usage: sh supertestscript.sh [ [ 0-3 ] | [raxml options] ]"
  echo " options: "
  echo "        [0]  tiny size test"
  echo "        [1]  small size test"
  echo "        [2]  medium size test"
  echo "        [3]  large size test"
  echo "        [raxml options]  specific test"
  echo "\n"
  exit 0
fi

if [ $# -eq 1 ] ; then
  if [ $1 -ne 0 ] && [ $1 -ne 1 ] && [ $1 -ne 2 ] && [ $1 -ne 3 ] ; then
    echo "\n"
    echo " usage: sh supertestscript.sh [ [ 0-3 ] | [raxml options] ]"
    echo " options: "
    echo "        [0]  tiny size test"
    echo "        [1]  small size test"
    echo "        [2]  medium size test"
    echo "        [3]  large size test"
    echo "        [raxml options]  specific test"
    echo "\n"
    exit 1
  fi
  if [ $1 -eq 0 ] ; then
    DESC="tiny"
  fi
  if [ $1 -eq 1 ] ; then
    DESC="small"
  fi
  if [ $1 -eq 2 ] ; then
    DESC="medium"
  fi
  if [ $1 -eq 3 ] ; then
    DESC="large"
  fi

  #SIMPLE=""                # this will be ignored in the loop
  #BL_PARTITION="-M"
  #SIMPLE_GAPPY="-S"
  #BL_PARTITION_GAPPY="-M -S" # !! this will be taken as 2
  #SIMPLE_RF_CONV="-D"

  echo "Starting recom custom test `date`, Errors" > $ERRLOGFILE 
  echo "Starting recom custom test `date`, Log" > $LOGFILE 

  for VERSION in SSE3_GCC SSE3_PTHREADS_GCC 
  do 
    echo $VERSION
    for MODEL in PSR GAMMA 
    do
      echo $MODEL
      for PARTITION_LABEL in singlegene.binary binary 
      do
        for ALPHABET in dna aa 
        do
          echo $ALPHABET
          echo "   "
          TREE="${DATADIR}/${DESC}.startingTree.${ALPHABET}.tree"
          ALIGNMENT="${DATADIR}/${DESC}.${ALPHABET}.${PARTITION_LABEL}"
          BASE_OPTIONS="-m ${MODEL} -s ${ALIGNMENT} -t ${TREE}"
          echo ${BASE_OPTIONS} 

          #Normal
          OPTIONS="${BASE_OPTIONS}"
          echo "$VERSION with $OPTIONS" | tee -a $ERRLOGFILE $LOGFILE
          (run_${VERSION} 2>> $ERRLOGFILE) >> $LOGFILE

          # RF convergence
          OPTIONS="${BASE_OPTIONS} -D"
          echo "$VERSION with $OPTIONS" | tee -a $ERRLOGFILE $LOGFILE
          (run_${VERSION} 2>> $ERRLOGFILE) >> $LOGFILE

        done
      done
    done
  done
else
  OPTIONS="$@"
  for VERSION in SSE3_GCC SSE3_PTHREADS_GCC 
  do 
    run_${VERSION}
  done
fi

echo "\n"
echo " recom test script finished successfully. Check file ${LOGFILE} to ensure that everything is ok!"
echo "\n"

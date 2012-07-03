#! /bin/sh
# This is the working directory you want to run the test from
WORKINGDIR=`pwd`
DATADIR=`pwd`/../testdata
# Set the number of threads you want to use
NUMPTHREADS=4
# NO NEED TO EDIT BEYOND THIS POINT
AVX_GCC="testscript.AVX.gcc"
AVX_PTHREADS_GCC="testscript.PTHREADS.AVX.gcc"
SSE3_GCC="testscript.SSE3.gcc"
SSE3_PTHREADS_GCC="testscript.SSE3.PTHREADS.gcc"

# KEEP A LOG OF STDOUT OF RAxML
LOGFILE=RAxML_log.txt
ERRLOGFILE=RAxML_ERR_log.txt

run_AVX_GCC()
{
  	sh ${WORKINGDIR}/${AVX_GCC} ${OPTIONS} 
  	RETVAL=$?
  	if [ $RETVAL -ne 0 ] ; then
  		echo " Error: testscript.AVX.gcc failed. Check the messages above ...";
		echo "\n"
		echo " Supertestscript did not finish successfully."
		echo "\n"
        	exit 1
  	fi
  	echo "\n"
}

run_AVX_PTHREADS_GCC()
{
	sh ${WORKINGDIR}/${AVX_PTHREADS_GCC} ${OPTIONS} -T ${NUMPTHREADS}
	RETVAL=$?
	if [ $RETVAL -ne 0 ] ; then
        	echo " Error: testscript.PTHREADS.AVX.gcc failed. Check the messages above ...";
		echo "\n"
		echo " Supertestscript did not finish successfully."
		echo "\n"
        	exit 1
	fi
	echo "\n"
}

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

if [ ! -f ${WORKINGDIR}/${AVX_GCC} ] ; then
        echo " Error: file ${AVX_GCC} is missing from ${WOKRINGDIR}, exiting ...";
        exit 1
fi
if [ ! -f ${WORKINGDIR}/${AVX_PTHREADS_GCC} ] ; then
        echo " Error: file ${AVX_PTHREADS_GCC} is missing from ${WORKINGDIR}, exiting ...";
        exit 1
fi
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
		TEST_DNA_PARTITIONED="${DATADIR}/tiny.dna.binary"
		TEST_DNA_SINGLE="${DATADIR}/tiny.dna.singlegene.binary"
		TEST_AA_PARTITIONED="${DATADIR}/tiny.aa.binary"
		TEST_AA_SINGLE="${DATADIR}/tiny.aa.singlegene.binary"
        fi
	if [ $1 -eq 1 ] ; then
		TEST_DNA_PARTITIONED="${DATADIR}/small.dna.binary"
		TEST_DNA_SINGLE="${DATADIR}/small.dna.singlegene.binary"
		TEST_AA_PARTITIONED="${DATADIR}/small.aa.binary"
		TEST_AA_SINGLE="${DATADIR}/small.aa.singlegene.binary"
        fi
	if [ $1 -eq 2 ] ; then
		TEST_DNA_PARTITIONED="${DATADIR}/medium.dna.binary"
		TEST_DNA_SINGLE="${DATADIR}/medium.dna.singlegene.binary"
		TEST_AA_PARTITIONED="${DATADIR}/medium.aa.binary"
		TEST_AA_SINGLE="${DATADIR}/medium.aa.singlegene.binary"
        fi
	if [ $1 -eq 3 ] ; then
		TEST_DNA_PARTITIONED="${DATADIR}/large.dna.binary"
		TEST_DNA_SINGLE="${DATADIR}/large.dna.singlegene.binary"
		TEST_AA_PARTITIONED="${DATADIR}/large.aa.binary"
		TEST_AA_SINGLE="${DATADIR}/large.aa.singlegene.binary"
        fi

	SIMPLE="" 
	BL_PARTITION="-M"
	SIMPLE_GAPPY="-S"
	BL_PARTITION_GAPPY="-M -S"
	SIMPLE_RF_CONV="-D"

  echo "Starting superscript `date`, Errors" > $ERRLOGFILE 
  echo "Starting superscript `date`, Log" > $LOGFILE 

	for VERSION in AVX_GCC AVX_PTHREADS_GCC SSE3_GCC SSE3_PTHREADS_GCC
  	do 
    		for MODEL in PSR GAMMA 
     		do
       			for DATA_PARTITIONED in ${TEST_DNA_PARTITIONED} ${TEST_AA_PARTITIONED}
         		do
            			for FLAGS in ${BL_PARTITION} ${BL_PARTITION_GAPPY} 
              			do
                			OPTIONS="-m ${MODEL} -s ${DATA_PARTITIONED} ${FLAGS}"
                      echo "$VERSION with $OPTIONS" | tee -a $ERRLOGFILE $LOGFILE
                      (run_${VERSION} 2>> $ERRLOGFILE) >> $LOGFILE
              			done
         		done
       			for DATA_SINGLE in ${TEST_DNA_SINGLE} ${TEST_AA_SINGLE}
         		do
            			for FLAGS in ${SIMPLE} ${SIMPLE_GAPPY} ${SIMPLE_RF_CONV} 
              			do
                			OPTIONS="-m ${MODEL} -s ${DATA_SINGLE} ${FLAGS}"
                      echo "$VERSION with $OPTIONS" | tee -a $ERRLOGFILE $LOGFILE
                      (run_${VERSION} 2>> $ERRLOGFILE) >> $LOGFILE
              			done
         		done
     		done
  	done
else
	OPTIONS="$@"
	for VERSION in AVX_GCC AVX_PTHREADS_GCC SSE3_GCC SSE3_PTHREADS_GCC
  	do 
		run_${VERSION}
        done
fi

echo "\n"
echo " Supertestscript finished successfully."
echo "\n"

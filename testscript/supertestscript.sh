#! /bin/sh
# This is the working directory you want to run the test from
WORKINGDIR="/home/pississn/Library/local/testscript"
DATADIR="/home/pississn/Library/local/testdata"
# NO NEED TO EDIT BEYOND THIS POINT
AVX_GCC="testscript.AVX.gcc"
AVX_PTHREADS_GCC="testscript.PTHREADS.AVX.gcc"
SSE3_GCC="testscript.SSE3.gcc"
SSE3_PTHREADS_GCC="testscript.SSE3.PTHREADS.gcc"

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
	sh ${WORKINGDIR}/${AVX_PTHREADS_GCC} ${OPTIONS}
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
	sh ${WORKINGDIR}/${SSE3_PTHREADS_GCC} ${OPTIONS}
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
	echo "        [0]  supertiny size test"
	echo "        [1]  tiny size test"
	echo "        [2]  medium size test"
	echo "        [3]  huge size test"
	echo "        [raxml options]  specific test"
	echo "\n"
	exit 0
fi

if [ $# -eq 1 ] ; then
        if [ $1 -ne 0 ] && [ $1 -ne 1 ] && [ $1 -ne 2 ] && [ $1 -ne 3 ] ; then
		echo "\n"
		echo " usage: sh supertestscript.sh [ [ 0-3 ] | [raxml options] ]"
		echo " options: "
		echo "        [0]  supertiny size test"
		echo "        [1]  tiny size test"
		echo "        [2]  medium size test"
		echo "        [3]  huge size test"
		echo "        [raxml options]  specific test"
		echo "\n"
		exit 1
        fi
	if [ $1 -eq 0 ] ; then
		TEST1="${DATADIR}/7.dna.binary"
		TEST2="${DATADIR}/7.dna.singlegene.binary"
		TEST3="${DATADIR}/7.aa.binary"
		TEST4="${DATADIR}/7.aa.singlegene.binary"
        fi
	if [ $1 -eq 1 ] ; then
		TEST1="${DATADIR}/20.dna.binary"
		TEST2="${DATADIR}/20.dna.singlegene.binary"
		TEST3="${DATADIR}/20.aa.binary"
		TEST4="${DATADIR}/20.aa.singlegene.binary"
        fi
	if [ $1 -eq 2 ] ; then
		TEST1="${DATADIR}/5K.dna.binary"
		TEST2="${DATADIR}/5K.dna.singlegene.binary"
		TEST3="${DATADIR}/5K.aa.binary"
		TEST4="${DATADIR}/5K.aa.singlegene.binary"
        fi
	if [ $1 -eq 3 ] ; then
		TEST1="${DATADIR}/40K.dna.binary"
		TEST2="${DATADIR}/40K.dna.singlegene.binary"
		TEST3="${DATADIR}/40K.aa.binary"
		TEST4="${DATADIR}/40K.aa.singlegene.binary"
        fi

	SIMPLE="" 
	BL_PARTITION="-M"
	SIMPLE_GAPPY="-S"
	BL_PARTITION_GAPPY="-M -S"
	SIMPLE_RF_CONV="-D"

	for VERSION in AVX_GCC AVX_PTHREADS_GCC SSE3_GCC SSE3_PTHREADS_GCC
  	do 
    		# DNA MODEL loop
    		for MODEL in PSR GAMMA 
     		do
       			for DATA_PARTITIONED in ${TEST2} 
         		do
            			for FLAGS in ${BL_PARTITION} ${BL_PARTITION_GAPPY} 
              			do
                			OPTIONS="-m ${MODEL} -s ${DATA_PARTITIONED} ${FLAGS}"
					run_${VERSION}
              			done
         		done
       			for DATA_SINGLE in ${TEST1} 
         		do
            			for FLAGS in ${SIMPLE} ${SIMPLE_GAPPY} ${SIMPLE_RF_CONV} 
              			do
                			OPTIONS="-m ${MODEL} -s ${DATA_SINGLE} ${FLAGS}"
					run_${VERSION}
              			done
         		done
     		done

    		# AA MODEL loop
#    		for MODEL in PROTPSRAUTO PROTPSRGAMMA PROTPSRAUTOF PROTGAMMAAUTOF
#     		do
#       			for DATA_PARTITIONED in ${TEST4} 
#         		do
#            			for FLAGS in ${BL_PARTITION} ${BL_PARTITION_GAPPY} 
#              			do
#                			OPTIONS="-m ${MODEL} -s ${DATA_PARTITIONED} ${FLAGS}"
#					run_${VERSION}
#              			done
#         		done
#       			for DATA_SINGLE in ${TEST3} 
#         		do
#            			for FLAGS in ${SIMPLE} ${SIMPLE_GAPPY} ${SIMPLE_RF_CONV} 
#              			do
#                			OPTIONS="-m ${MODEL} -s ${DATA_SINGLE} ${FLAGS}"
#					run_${VERSION}
#              			done
#         		done
#     		done
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

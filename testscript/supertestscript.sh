#! /bin/sh
# This is the working directory you want to run the test from
WORKINGDIR="/home/pississn/Library/local/testscript"

if [ ! -f ${WORKINGDIR}/testscript.AVX.gcc ] ; then
        echo " Error: file testscript.AVX.gcc is missing from ${WOKRINGDIR}, exiting ...";
        exit 0
fi
if [ ! -f ${WORKINGDIR}/testscript.PTHREADS.AVX.gcc ] ; then
        echo " Error: file testscript.PTHREADS.AVX.gcc is missing from ${WORKINGDIR}, exiting ...";
        exit 0
fi
if [ ! -f ${WORKINGDIR}/testscript.SSE3.gcc ] ; then
        echo " Error: file testscript.SSE3.gcc is missing from ${WORKINGDIR}, exiting ...";
        exit 0
fi
if [ ! -f ${WORKINGDIR}/testscript.SSE3.PTHREADS.gcc ] ; then
        echo " Error: file testscript.SSE3.PTHREADS.gcc is missing from ${WORKINGDIR}, exiting ...";
        exit 0
fi

if [ $# -eq 0 ] ; then
	echo "\n"
	echo " usage: sh supertestscript.sh [raxml options]"
	echo "\n"
	exit 1
fi

sh ${WORKINGDIR}/testscript.AVX.gcc $@
RETVAL=$?
if [ $RETVAL -ne 0 ] ; then
        echo " Error: testscript.AVX.gcc failed. Check the messages above ...";
	echo "\n"
	echo " Supertestscript did not finish successfully."
	echo "\n"
        exit 1
fi
echo "\n"
sh ${WORKINGDIR}/testscript.PTHREADS.AVX.gcc $@
RETVAL=$?
if [ $RETVAL -ne 0 ] ; then
        echo " Error: testscript.PTHREADS.AVX.gcc failed. Check the messages above ...";
	echo "\n"
	echo " Supertestscript did not finish successfully."
	echo "\n"
        exit 1
fi
echo "\n"
sh ${WORKINGDIR}/testscript.SSE3.gcc $@
RETVAL=$?
if [ $RETVAL -ne 0 ] ; then
        echo " Error: testscript.SSE3.gcc failed. Check the messages above ...";
	echo "\n"
	echo " Supertestscript did not finish successfully."
	echo "\n"
        exit 1
fi
echo "\n"
sh ${WORKINGDIR}/testscript.SSE3.PTHREADS.gcc $@
RETVAL=$?
if [ $RETVAL -ne 0 ] ; then
        echo " Error: testscript.SSE3.PTHREADS.gcc failed. Check the messages above ...";
	echo "\n"
	echo " Supertestscript did not finish successfully."
	echo "\n"
        exit 1
fi
echo "\n"

echo "\n"
echo " Supertestscript finished successfully."
echo "\n"

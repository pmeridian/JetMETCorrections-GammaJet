#!/bin/sh

#set -x
###FIXME CHECK THAT 2 STRINGS ARE GIVEN IN INPUT
sleeptime=1200

sleepAndRenew(){
    echo "Sleeping $sleeptime"
    sleep $sleeptime
}

usage(){
    echo "$0 [friendly dataset name] [output dir in EOS]"
    exit
}

if [ $# -lt 2 ] 
then
    usage
fi

if [ ! $CRABDIR ] ; then
   echo "Please set CRAB environment"
   exit
fi

if [ ! $CRABDIR ] ; then
   echo "Please set CRAB environment"
   exit
fi
#source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh
#eval `scramv1 runtime -sh`
#source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh

STARTDIR=`pwd`
CASTORDIR=$2
#BLACKLIST=`cat ~/physEGAMMA/dashboard/blacklist.txt` 
while [ "X" = "X" ]; do
    cd ${STARTDIR}
    for run in `find  $1 -maxdepth 1 -type d | sed 1d | grep -v ".old" | sed -e "s%\./%%g"`;  do
	pwd
	run="${run%\\n}"
	cd ${STARTDIR}
	task=${run}
#	cd $task
	echo "[`date`]: ++++++++++ Checking task $task ++++++ "
#	ls $run
#    rm -rf ${task}/res/lumiSummary.json
	
	#crabdir=`find . -type d -iname "crab_*" | tail -1`
	crabdir=$run
	if [ -f ${crabdir}/res/lumiSummary.json ]; then
	    echo ">>>> Task ${task} already cleared. Moving on"
	    continue
	fi
	
    #CHECKING STATUS. HAVE TO DO IT TWICE SOMETIMES CRAB SERVER IS CRAZY!!
#	echo "crab -status -c $crabdir > ${task}_status.log"
	crab -status -c $crabdir > ${task}_status.log
#	sleep 2
#	crab -status -c $crabdir > ${task}_status.log 2>&1
	
	nTotalJobs=`cat ${task}_status.log | grep "Total Jobs" | awk '{print $2}'`
	jobsSubmittingCreated=`cat ${task}_status.log | egrep "(Submit|Created|Running|RUN|Ready)" `
	if [ "${jobsSubmittingCreated}AAA" != "AAA" ]; then
	    echo ">>>> ${task} still in submission/created/running/ready state. Better to wait or relaunch it.\n>>>> Do you want to relaunch it now (y/n)? 10 seconds to decide"
#	    read -t 10 KEYINPUT
#	    if [ "$KEYINPUT" == "y" ]; then
#		eos rm -r ${CASTORDIR}/${task}; mv -f $task $task.old; crab -create -submit --GRID.ce_black_list='ce02.lcg.cscs.ch,ce11.lcg.cscs.ch' --GRID.se_black_list=${BLACKLIST}  -cfg ${task}.cfg > ${task}_resubmit.log 2>&1 &
#	    fi
	    continue
	fi

	crab -getoutput ${jobsCrashed} -c $crabdir > ${task}_getoutput.log
#	jobsDoneOK=`cat ${task}_status.log | egrep -i "(Done)" | awk '{printf "%d,",$1}' | sed -e 's%,$%%g'`
	jobsDoneOK=`cat ${task}_status.log | grep -i -A 2 "Wrapper Exit Code" | grep -A 2 "Wrapper Exit Code : 0" | grep -A2 ">>>>>>>>>" | grep "List of jobs:" | sed -e "s%.*List of jobs: %%g" | awk '{printf "%s,",$1}' | sed -e 's%,$%%g' | awk -f ${STARTDIR}/expandJobsList.awk | sed -e 's%,$%%g'`
	if [ "${jobsDoneOK}AAA" != "AAA" ]; then
	    echo ">>>> Jobs ${jobsDoneOK} are OK for task $task"
	fi

	jobsDoneDup=`cat ${task}_status.log | grep -i -A 2 "Wrapper Exit Code" | grep -A 2 "Wrapper Exit Code : 60303" | grep -A2 ">>>>>>>>>" | grep "List of jobs:" | sed -e "s%.*List of jobs: %%g" | awk '{printf "%s,",$1}' | sed -e 's%,$%%g' | awk -f ${STARTDIR}/expandJobsList.awk | sed -e 's%,$%%g' ` 
	if [ "${jobsDoneDup}AAA" != "AAA" ]; then
	    echo ">>>> Jobs ${jobsDoneDup} are in status 60303 (FileAlready present in the output). In a future version can will be cleaned and relaunched"
	fi
	
	jobsCrashed=`cat ${task}_status.log | grep -i -A 2 "Wrapper Exit Code" | grep -v "Wrapper Exit Code : 0" | grep -v "Wrapper Exit Code : 60303" | grep -A2 ">>>>>>>>>" | grep "List of jobs:" | sed -e "s%.*List of jobs: %%g" | awk '{printf "%s,",$1}' | sed -e 's%,$%%g' | awk -f ${STARTDIR}/expandJobsList.awk | sed -e 's%,$%%g' ` 
	if [ "${jobsCrashed}AAA" != "AAA" ]; then
	    echo ">>>> Relaunching crashed jobs for task $task: ${jobsCrashed}"
	     crab -forceResubmit ${jobsCrashed} -c $crabdir > ${task}_resubmit_crashed.log&
	fi
	
	jobsAborted=`cat ${task}_status.log | grep -i -A 2 "Aborted" |  grep -A2 ">>>>>>>>>" | grep "List of jobs:" | sed -e "s%.*List of jobs: %%g" | awk '{printf "%s,",$1}' | sed -e 's%,$%%g' | awk -f ${STARTDIR}/expandJobsList.awk | sed -e 's%,$%%g'` 
	if [ "${jobsAborted}AAA" != "AAA" ]; then
	    echo ">>>> Relaunching aborted jobs for task $task: jobs ${jobsAborted}"
	    crab -forceResubmit ${jobsAborted}  -c $crabdir > ${task}_resubmit_aborted.log&
	fi
	
	
    ##FINAL CHECKS. IF EVERYTHING OK THEN CLEARING OUTPUT. CONSIDERING GOOD FOR THE MOMENT ALSO 60303
	jobsDone=${jobsDoneOK},${jobsDoneDup},
	outDir=`echo ${task} | awk -F '/' '{print $2}'`
	if [ "${jobsDone}AAA" != ",AAA," ]; then
	    echo ">>>> Checking output of jobs for task $task"
	    FILENOTOK=""
	    for job in `echo ${jobsDone} | sed -e 's%,% %g'`
	      do
	      FILEOK=`/afs/cern.ch/project/eos/installation/0.1.0-22d/bin/eos.select ls ${CASTORDIR}/${outdir} | grep  output_${job} | awk '{print $9}'`
#	      if [ "${FILEOK}AAA" == "AAA" ]; then
#	       	  FILENOTOK=`echo "${FILENOTOK}${job},"`
#	      fi
	    done
	    if [ "${FILENOTOK}AAA" != "AAA" ]; then
		echo ">>>> Jobs ${FILENOTOK} do not have their output. Relaunching them"
		jobsCrashed=`echo ${FILENOTOK} | sed -e 's%,$%%g'` 
		crab -forceResubmit ${jobsCrashed} -c $crabdir > ${task}_resubmit_crashed.log&
	    else
		nJobsDoneOK=`echo $jobsDoneOK | awk -F "," '{print NF}'`
		nJobsDoneDup=`echo $jobsDoneDup | awk -F "," '{print NF}'`
		(( nJobsDone = nJobsDoneOK + nJobsDoneDup ))
		if [ "${jobsCrashed}AAA" == "AAA"  -a "${jobsAborted}AAA" == "AAA" -a ${nTotalJobs} == ${nJobsDone} ]; then
		    echo ">>>> Everything OK. Clearing task ${task}"
		    crab -report -c $crabdir > ${task}_report.log&
		    tar cvzf ${task}_getoutput.tar.gz $crabdir/res/CMSSW*; rm -rf $crabdir/res/CMSSW*; 
#mv /tmp/meridian/${task}_getoutput.tar.gz ./
		else
		    echo ">>>> Still not all jobs done/retrieved/ok for ${task}. Try later"
		fi
	    
	    fi
	fi
    done
    sleepAndRenew
done
#!/bin/sh

subDir=`pwd`

outDir=/eos/uscms/store/user/meridian/MC/53xv4/
run=0
check=0

while getopts d:rc flag; do
  case $flag in
    r)
      echo "Submit is on";
      run=1
      ;;
    c)
      echo "Check is on";
      check=1
      ;;
    d)
      echo "outdir is $OPTARG";
      outDir=$OPTARG
      ;;
    ?)
      exit;
      ;;
  esac
done
shift $(( OPTIND - 1 ));

if [ "$run" == "1" ]; then
    for dir in `find $1 -maxdepth 1 -type d`; do
	echo "Submitting ${dir}"
#	cd ${dir}/share/
#	mkdir -p tt
#	cd tt/
#	echo "Patching  ${dir}/share/default.tgz"
#	tar xvzf ../default.tgz > /dev/null
#	cat cmscp.py | sed -e "s%3600%7200%g" > cmscp.py.tmp > /dev/null
#	mv -f cmscp.py.tmp cmscp.py > /dev/null
#	tar cvzf ../default.tgz * > /dev/null
#	cd  ..
#	rm -rf tt/
#	cd ${subDir}
#	ls -lah ${dir}/share/default.tgz
	./submitTask.py -t ${dir} -s 5 -n 500 > ${dir}_submit.log 2>&1 
    done
fi

if [ "$check" == "1" ]; then
    echo "Launching checkJobs daemon for tasks in ${1} writing output in ${outDir}"
    ./checkAndRelaunchJobs.sh ${1} ${outDir} > checkJobs_${1}.log 2>&1 &
fi

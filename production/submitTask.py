#!/usr/bin/env python

from optparse import OptionParser, make_option
import sys
import os
import math
import subprocess as sub

def getNumberOfJobs(dir):
#    command='cat '+dir+'/log/crab.log | grep "Total jobs" | awk -F "Total jobs " \'{print $2}\' | head -n 1'
    command='cat '+dir+'/log/crab.log | grep "Total of" | awk -F "Total of " \'{print $2}\' | awk \'{print $1}\' | head -n 1'
    print command
    p = sub.Popen(command,shell=True,stdout=sub.PIPE,stderr=sub.PIPE)
    njobs=0
    for line in iter(p.stdout.readline,''):
        print line.rstrip()
        njobs=line.rstrip()
    return njobs

def submitJobs(dir,startJob,endJob,sleepTime,resubmit):
    command='crab '
    if (resubmit):
        command+='-forceResubmit '
    else:
        command+='-submit '
    command+=str(startJob)+'-'+str(endJob)+' -c '+dir
    sub.call(command,shell=True,stderr=sub.STDOUT)
    print "-_- Sleeping "+str(sleepTime)+"s -_-" 
    sub.call('sleep '+str(sleepTime),shell=True,stderr=sub.STDOUT)


def submitJobsWithList(dir,list,sleepTime,resubmit):
    command='crab '
    if (resubmit):
        command+='-forceResubmit '
    else:
        command+='-submit '
    command+=list+' -c '+dir
    sub.call(command,shell=True,stderr=sub.STDOUT)
    print command
    print "-_- Sleeping "+str(sleepTime)+"s -_-" 
    sub.call('sleep '+str(sleepTime),shell=True,stderr=sub.STDOUT)

def main(options,args):
    nJobs=0
    if (options.jobList == ""):
        nJobs=getNumberOfJobs(options.taskDir)
        print "The task "+options.taskDir+ " has "+str(nJobs)+" jobs" 
        for i in range(0,int(math.floor((int(nJobs)-int(options.firstJob)-1)/options.nJobs)+1)):
            print "+++ submitting jobs "+str(int(i*options.nJobs+1))+"-"+str(int((i+1)*options.nJobs)) 
            submitJobs(options.taskDir,int(i*options.nJobs+options.firstJob),int((i+1)*options.nJobs+options.firstJob-1),options.sleepTime,options.resubmit)
    else:
        jobs=options.jobList.split(',')
        nJobs=len(jobs)
        print "Launching for task "+options.taskDir+ " n: "+str(nJobs)+" jobs" 
        #        print jobs
        for i in range(0,int(math.floor((int(nJobs)-1)/options.nJobs)+1)):
            submitJobList=''
            for job in range(int(i*options.nJobs),min(int((i+1)*options.nJobs),len(jobs))):
                submitJobList+=str(jobs[job])+","
            submitJobList=submitJobList[:-1]
            print "+++  submitting jobs "+str(submitJobList)
            submitJobsWithList(options.taskDir,submitJobList,options.sleepTime,options.resubmit)
        
if __name__ == "__main__":
    parser = OptionParser(option_list=[
        make_option("-t", "--taskDir",
                    action="store", type="string", dest="taskDir",
                    default="",
                    help="", metavar=""
                    ),
        make_option("-s", "--sleepTime",
                    action="store", type="int", dest="sleepTime",
                    default=3600,
                    help="", metavar=""
                    ),
        make_option("-n", "--nJobs",
                    action="store", type="int", dest="nJobs",
                    default=500,
                    help="", metavar=""
                    ),
        make_option("-f", "--firstJob",
                    action="store", type="int", dest="firstJob",
                    default=1,
                    help="", metavar=""
                    ),
        make_option("-l", "--list",
                    action="store", type="string", dest="jobList",
                    default="",
                    help="", metavar=""
                    ),
        make_option("-r", "--resubmit",
                    action="store_true",  dest="resubmit",
                    help="", metavar=""
                    ),
        ])
    
    (options, args) = parser.parse_args()
    
    main( options, args) 


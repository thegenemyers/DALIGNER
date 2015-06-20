from pypeflow.common import *
from pypeflow.data import PypeLocalFile, makePypeLocalFile, fn
from pypeflow.task import PypeTask, PypeThreadTaskBase, PypeTaskBase
from pypeflow.controller import PypeWorkflow, PypeThreadWorkflow
import os
import uuid
import sys


def run_script(job_data, job_type = "SGE" ):
    if job_type == "SGE":
        job_name = job_data["job_name"]
        cwd = job_data["cwd"]
        sge_option = job_data["sge_option"]
        script_fn = job_data["script_fn"]
        sge_cmd="qsub -N {job_name} {sge_option} -o {cwd}/sge_log -j y\
                 -S /bin/bash {script}".format(job_name=job_name,
                                               cwd=os.getcwd(),
                                               sge_option=sge_option,
                                               script=script_fn)

        #print sge_cmd
        os.system( sge_cmd )
        os.system( "sleep 1")
    elif job_type == "local":
        os.system( "bash %s" % job_data["script_fn"] )

def wait_for_file(filename, task = None, job_name = ""):
    while 1:
        time.sleep(60)
        if os.path.exists(filename):
            break

        if task != None:
            if task.shutdown_event != None and task.shutdown_event.is_set():
                os.system("qdel %s" % job_name)
                break

def run_daligner(self):
    daligner_cmd = self.parameters["daligner_cmd"]
    job_id = self.parameters["job_id"]
    cwd = self.parameters["cwd"]
    script_dir = os.path.join( cwd )
    script_fn =  os.path.join( script_dir , "rj_%05d.sh" % (job_id))
    log_path = os.path.join( script_dir, "rj_%05d.log" % (job_id))
    script = []
    script.append( "export PATH=~/task2014/dazzler/DALIGNER/:$PATH" )
    script.append( "cd %s" % cwd )
    script.append( "/usr/bin/time "+ daligner_cmd + ( " >& %s " % log_path ) + ( " && touch %s" % fn( self.job_done ) ) )

    with open(script_fn,"w") as script_file:
        script_file.write("\n".join(script))

    job_name = self.URL.split("/")[-1]
    job_name += "-"+str(uuid.uuid1())[:8]
    job_data = {"job_name": job_name,
                "cwd": cwd,
                "sge_option": " -pe smp 6 -q huasm ",
                "script_fn": script_fn }
    run_script(job_data, job_type = "SGE")
    wait_for_file( fn( self.job_done ), task=self, job_name=job_name )

if __name__ == "__main__":
    prefix = sys.argv[1]
    concurrent_jobs = 64
    PypeThreadWorkflow.setNumThreadAllowed(concurrent_jobs, concurrent_jobs)
    wf = PypeThreadWorkflow()

    job_id = 0
    db_file = makePypeLocalFile(os.path.abspath( "./%s.db" % prefix ))
    with open("run_jobs.sh") as f :
        for l in f :
            l = l.strip().split()
            if l[0] == "daligner":
                try:
                    os.makedirs("./job_%05d" % job_id)
                except OSError:
                    pass
                os.system("cd ./job_%05d;ln -s ../.%s.bps .; ln -s ../.%s.idx .; ln -s ../%s.db ." % (job_id, prefix, prefix, prefix) )
                job_done = makePypeLocalFile(os.path.abspath( "./job_%05d/job_%05d_done" % (job_id,job_id)  ))
                parameters =  {"daligner_cmd": " ".join(l),
                               "cwd": os.path.join(os.getcwd(), "job_%05d" % job_id),
                               "job_id": job_id}
                make_daligner_task = PypeTask( inputs = {"db_file": db_file},
                                               outputs = {"job_done": job_done},
                                               parameters = parameters,
                                               TaskType = PypeThreadTaskBase,
                                               URL = "task://localhost/mtask_%05d" % job_id )
                daligner_task = make_daligner_task ( run_daligner )
                wf.addTask(daligner_task)
                job_id += 1
                print job_id
    wf.refreshTargets(updateFreq = 45) #all

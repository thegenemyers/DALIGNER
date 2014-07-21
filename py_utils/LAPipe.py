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
        time.sleep(30)
        if os.path.exists(filename):
            break

        if task != None:
            if task.shutdown_event != None and task.shutdown_event.is_set():
                os.system("qdel %s" % job_name)
                break

def run_p_task(self):
    p_script_fn = self.parameters["p_file"]
    job_id = self.parameters["job_id"]
    cwd = self.parameters["cwd"]
    script_dir = os.path.join( cwd )
    script_fn =  os.path.join( script_dir , "rp_%05d.sh" % (job_id))
    log_path = os.path.join( script_dir, "rp_%05d.log" % (job_id))
    script = []
    script.append( "export PATH=~/task2014/dazzler/DALIGNER/:$PATH" )
    script.append( "cd %s" % cwd )
    script.append( ("/usr/bin/time bash %s " % p_script_fn)  + ( " >& %s " % log_path ) + ( " && touch %s" % fn( self.job_done ) ) )

    with open(script_fn,"w") as script_file:
        script_file.write("\n".join(script))



    job_name = self.URL.split("/")[-1]
    job_name += "-"+str(uuid.uuid1())[:8]
    job_data = {"job_name": job_name,
                "cwd": cwd,
                "sge_option": " -pe smp 2 -q huasm ",
                "script_fn": script_fn }
    run_script(job_data, job_type = "SGE")
    wait_for_file( fn( self.job_done ), task=self, job_name=job_name )

def run_consensus_task(self):
    job_id = self.parameters["job_id"]
    cwd = self.parameters["cwd"]
    script_dir = os.path.join( cwd )
    script_fn =  os.path.join( script_dir , "cp_%05d.sh" % (job_id))
    log_path = os.path.join( script_dir, "cp_%05d.log" % (job_id))

    with open( os.path.join(cwd, "c_%05d.sh" % job_id), "w") as p_script:
        print >> p_script, ". /mnt/secondary/Share/HBAR_03202013/bin/activate"
        print >> p_script, "cd .."
        print >> p_script, """./LA4Falcon -o -f:%s las_files/%s.%d.las | """ % (prefix, prefix, job_id),
        print >> p_script, """ falcon_sense.py --trim --output_multi --min_idt 0.70 --min_cov 4 --local_match_count_threshold 3 --max_n_read 800 --n_core 8 > %s""" % fn(self.out_file)

    script = []
    script.append( "cd %s" % cwd )
    script.append( ("/usr/bin/time bash c_%05d.sh " % job_id )  + ( " >& %s " % log_path ) + ( " && touch c_%05d_done" % job_id  ) )

    with open(script_fn,"w") as script_file:
        script_file.write("\n".join(script))

    job_name = self.URL.split("/")[-1]
    job_name += "-"+str(uuid.uuid1())[:8]
    job_data = {"job_name": job_name,
                "cwd": cwd,
                "sge_option": " -pe smp 6 -q huasm ",
                "script_fn": script_fn }
    run_script(job_data, job_type = "SGE")
    wait_for_file( os.path.join(cwd,"c_%05d_done" % job_id) , task=self, job_name=job_name )


if __name__ == "__main__":

    prefix = sys.argv[1]

    concurrent_jobs = 16
    PypeThreadWorkflow.setNumThreadAllowed(concurrent_jobs, concurrent_jobs)
    wf = PypeThreadWorkflow()

    mjob_data = {}

    with open("run_jobs.sh") as f:
        for l in f:
            l = l.strip().split()
            if l[0] not in ( "LAsort", "LAmerge" ):
                continue
            if l[0] == "LAsort":
                p_id = int( l[2].split(".")[1] )
                mjob_data.setdefault( p_id, [] )
                mjob_data[p_id].append(  " ".join(l) )
            if l[0] == "LAmerge":
                l2 = l[2].split(".")
                if l2[1] == "L2":
                    p_id = int(  l[2].split(".")[2] )
                    mjob_data.setdefault( p_id, [] )
                    mjob_data[p_id].append(  " ".join(l) )
                else:
                    p_id = int( l[2].split(".")[1] )
                    mjob_data.setdefault( p_id, [] )
                    mjob_data[p_id].append(  " ".join(l) )

    db_file = makePypeLocalFile(os.path.abspath( "./%s.db" % prefix ))

    for p_id in mjob_data:
        s_data = mjob_data[p_id]

        try:
            os.makedirs("./p_%05d" % p_id)
            os.makedirs("./p_%05d/sge_log" % p_id)
        except OSError:
            pass
        try:
            os.makedirs("./preads")
        except OSError:
            pass
        try:
            os.makedirs("./las_files")
        except OSError:
            pass
        with open("./p_%05d/p_%05d.sh" % (p_id, p_id), "w") as p_script:
            print >> p_script, """for f in `find .. -wholename "*job*/%s.%d.%s.*.*.las"`; do ln -sf $f .; done""" % (prefix, p_id, prefix)
            for l in s_data:
                print >> p_script, l
                print >> p_script, "mv %s.%d.las ../las_files" % (prefix, p_id)

        p_file = os.path.abspath( "./p_%05d/p_%05d.sh" % (p_id, p_id) )
        job_done = makePypeLocalFile(os.path.abspath( "./p_%05d/p_%05d_done" % (p_id,p_id)  ))
        parameters =  {"p_file": p_file,
                       "cwd": os.path.join(os.getcwd(), "p_%05d" % p_id),
                       "job_id": p_id}
        make_p_task = PypeTask( inputs = {"db_file": db_file},
                                       outputs = {"job_done": job_done},
                                       parameters = parameters,
                                       TaskType = PypeThreadTaskBase,
                                       URL = "task://localhost/ptask_%05d" % p_id )
        p_task = make_p_task ( run_p_task )

        wf.addTask(p_task)


        out_file = makePypeLocalFile(os.path.abspath( "./preads/out.%04d.fa" % p_id  ))
        parameters =  {"cwd": os.path.join(os.getcwd(), "preads" ),
                       "job_id": p_id}
        make_c_task = PypeTask( inputs = {"job_done": job_done},
                                outputs = {"out_file": out_file },
                                parameters = parameters,
                                TaskType = PypeThreadTaskBase,
                                URL = "task://localhost/ct_%05d" % p_id )

        c_task = make_c_task( run_consensus_task )
        wf.addTask(c_task)
        print p_id
    wf.refreshTargets(updateFreq = 15) #all

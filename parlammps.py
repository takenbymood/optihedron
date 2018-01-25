from multiprocessing import Pool, TimeoutError, Process
import time
import os
import subprocess
import plammps

from lammps import lammps
import signal

n = 1
nWorkers = 8



class Signal(Exception):
    """
    This exception is raise by the signal handler.
    """
    pass


class Timeout(Exception):
    """
    This exception is raised when the command exceeds the defined timeout
    duration and the command is killed.
    """
    def __init__(self, cmd, timeout):
        self.cmd = cmd
        self.timeout = timeout

    def __str__(self):
        return "Command '%s' timed out after %d second(s)." % \
               (self.cmd, self.timeout)


class Retcode(Exception):
    """
    This exception is raise when a command exits with a non-zero exit status.
    """
    def __init__(self, cmd, retcode, output=None):
        self.cmd = cmd
        self.retcode = retcode
        self.output = output

    def __str__(self):
        return "Command '%s' returned non-zero exit status %d" % \
               (self.cmd, self.returncode)


def alarm_handler(signum, frame):
    raise Signal


def execute(cmd, timeout=None):
    """
    Execute a command in the default shell. If a timeout is defined the command
    will be killed if the timeout is exceeded and an exception will be raised.
    Inputs:
        cmd     (str): Command to execute
        timeout (int): Command timeout in seconds
    Outputs:
        output (str): STDOUT/STDERR
    """
    # Define the timeout signal
    if timeout:
        signal.signal(signal.SIGALRM, alarm_handler)
        signal.alarm(timeout)

    try:
        # Execute the command and wait for the subprocess to terminate
        # STDERR is redirected to STDOUT
        phandle = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, preexec_fn=os.setsid)

        # Read the stdout/sterr buffers and retcode
        output, _ = phandle.communicate()
        retcode = phandle.poll()
    except Signal:
        # Kill the running process
        phandle.kill()
        raise Timeout(cmd=cmd, timeout=timeout)
    except:
        raise
    else:
        # Possible race condition where alarm isn't disabled in time
        signal.alarm(0)

    # Raise an exception if the command exited with non-zero exit status
    if retcode:
        raise Retcode(cmd, retcode, output=output)

    return output

def runSim(script,np,timeout):
    try:
        execute(['mpirun','-np',str(np),'python','./plammps.py','-s',str(script)],timeout)
    except:
        print('Process timed out')

def runSimSerial(script):
    try:
        plammps.startScript(script)
    except:
        print('Process crashed')

def runSims(scripts,np,timeout):
    processes = []
    for s in scripts:
        print(s)
        p = Process(target=runSim, args=(s,np,timeout,))
        p.start()
        processes.append(p)

    for p in processes:
        p.join()

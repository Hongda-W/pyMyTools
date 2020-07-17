"""
Print start, end, elapsed times
"""
import atexit
from time import time, strftime, localtime
from datetime import timedelta

def sec2str(elapsed=None):
    """ Pring duration in seconds to HH:MM:SS, current time will be returned if elapsed is None
    """
    if elapsed is None:
        return strftime("%Y-%m-%d %H:%M:%S", localtime())
    else:
        return str(timedelta(seconds=elapsed))

def log(s, elapsed=None):
    """ Write log information
    """
    line = '='*60
    print(line)
    print(sec2str(), '--', s)
    if elapsed is not None:
        print("Elapsed time:", elapsed)
    print(line)
    print()
    pass

def exitlog(start, func=''):
    """ Write log information at exit
    """
    end = time()
    elapsed = end-start
    log("End Program "+func, sec2str(elapsed))
    pass
    
def runtime(func=''):
    """ Get running time information for program, run time information will be written out upon normal interpreter termination
    """
    start = time()
    atexit.register(exitlog, start, func)
    log("Start Program "+func)
    return

def runtimef(func, *args, **kwargs):
    """ Get running time information for function
    """
    if not callable(func):
        raise ValueError("'%s' is not a function"%(func))
    start = time()
    log("Start Program "+func.__name__)
    func(*args, **kwargs)
    exitlog(start, func.__name__)
    return

""" Python I/O tools
"""

import numpy as np

def print2D(data, fmt=None):
    """ Print 2D array (list) to stdout, used fmt ("{:d} {:f}"...) to specify format
    """
    arr = np.asarray(data)
    for row in arr:
        if fmt is None:
            print(*row)
        else:
            print(fmt.format(*row))
    return

def write2D(data, outname=None, fmt=None, header=None):
    """ Write 2D array (list) to file if None print to stdout, used fmt ("{:d} {:f}"...) to specify format
    """
    if outname is None:
        print2D(data, fmt)
    arr = np.asarray(data)
    outf = open(outname, 'w+')
    if header is not None:
        outf.writelines(header.rstrip()+'\n')
    for row in arr:
        if fmt is not None:
            fmt = fmt.rstrip()+'\n'
            outf.writelines(fmt.format(*row))
        else:
            outf.writelines(" ".join(item for item in row)+'\n')
    return

""" Small system tools
"""

import pathlib
import shutil
import sys
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def mkdir(path, verbose=False):
    pathlib.Path(path).mkdir(parents=True, exist_ok=True)
    if verbose:
        print("Directory '{:s}' made".format(path))
    pass

def prompt_query(q, default=None):
    valid = {"yes": True, "y": True, "yep": True, "ay": True,
            "no": False, "n": False, "nay": False}
    if default is None:
        prompt = "[y/n]"
    elif default == "yes":
        prompt = "[Y/n]"
    elif default == "no":
        prompt = "[y/N]"
    else:
        raise ValueError("Invalid default answer: '{:s}'".format(default))
    while True:
        sys.stdout.write(q+prompt)
        choice = input().lower()
        if default is not None and choice=='':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no'"
            "(or 'y' or 'n')")

def rmdirF(path, prompt=True):
    """ Dangerous, rm -rf equivalent
    """
    if prompt:
        alarm = "ALARM! '{:s}' and its subdirectories will be removed!\n".format(path)
        sys.stdout.write(alarm)
        if prompt_query("Are you sure you want to continue?", default='no'):
            shutil.rmtree(path)
    pass
"""The Command-Line Interface (CLI) of sde
The CLI of sde can be accessed via ``python -m sde``.
:Example:
    Get help:
    
    .. code-block:: bash
        python -m sde -h
    
    Check version and authors:
    
    .. code-block:: bash
    
        python -m sde --version 
        python -m sde --author
**Pros**
    * Shallow learning curve
    
**Cons**
    * Fewer configurations 
    * No inspections of intermediate results
"""

import os
import sys
import argparse

from sde import __version__, __author__
from sde.sde_class import sde_class


parser = argparse.ArgumentParser(description="pMTnet Omni Document")

parser.add_argument("--version", action="version",
                    version=__version__, help="Display the version of the software")
parser.add_argument("--author", action="version", version=__author__,
                    help="Check the author list of the algorithm")

def main(cmdargs: argparse.Namespace):
    """The main method for sde

    Parameters:
    ----------
    cmdargs: argparse.Namespace
        The command line argments and flags 
    """

    sys.exit(0)


if __name__ == "__main__":
    cmdargs = parser.parse_args()
    main(cmdargs=cmdargs)

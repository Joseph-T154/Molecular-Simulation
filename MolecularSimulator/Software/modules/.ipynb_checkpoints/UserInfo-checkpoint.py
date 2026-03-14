#
# This script contains code that prints error, warning and info messages
# The weird text in each function is a color code
#
from enum import Enum

class ColorCodes(Enum):

    # Colors taken from some stackexchange page
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def printSuccess(msg):

    errStr = ColorCodes.OKBLUE.value + "OK: " + ColorCodes.ENDC.value + msg
    print(errStr)

def printWarning(msg):

    errStr = ColorCodes.WARNING.value + "Warning: " + ColorCodes.ENDC.value + msg
    print(errStr)
    
def printError(msg):

    errStr = ColorCodes.FAIL.value + "Error: " + ColorCodes.ENDC.value + msg
    print(errStr)
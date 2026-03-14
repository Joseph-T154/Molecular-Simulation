#
# This module contains all of the physical objects that can't easily be stored as lists
#

class LinearSpring:

    def __init__(self, k = [0], l = [0]):

        self.k = k
        self.l = l
        self.n = [-1,-1]
        self.currentState = 0
        

    def setIndices(self, a,b):

        self.n[0] = a
        self.n[1] = b
        
        
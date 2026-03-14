#
# This class reads in measurement data for analysis after the simulation has completed!
# This file just stores every measureable as a big list
# one element for each frame. Dead easy
#

import os
import numpy as np

class MDLJMeasurement:
    
    def __init__(self):
        
        self.simulationType = ""
        
        self.numFrames = 0
        self.time = None
        self.centroid = None
        self.kineticEnergy = None
        self.elasticEnergy = None
        self.LJEnergy = None
        
        self.numStages = 0
        self.stage = None
        
    def load(self, mFname):
        
        # Check filename exists
        if not os.path.exists(mFname):
            print("Cannot read, file does not exist :(")
            raise FileNotFoundError
        
        # Open and check it's the correct file
        fin = open(mFname, "r")

        line = fin.readline().strip()
        if line != "# MS LJGAS MEASUREMENT #":
            print("Cannot read, not the correct filetype :(")
            raise FileNotFoundError

        # Ok it should work from here. Any errors are probably Ben's fault...
        
        # Get simulation type
        self.simulationType = fin.readline().split(",")[1].strip()
            
        # Get headers
        headers = [bit.strip() for bit in fin.readline().split(",")]
        
        # Read header separator
        fin.readline()
        
        # Get number of frames before reading properly
        filePos = fin.tell()
        while True:
            line = fin.readline().strip()
            if line == "":
                break
            elif "##" in line:
                
                self.numStages += 1
                
                # Skip the stage marker
                for i in range(3):
                    fin.readline()
                
            self.numFrames += 1
            
        # Initialise memory structures
        self.time = np.zeros(self.numFrames, dtype=float)
        self.centroid = np.zeros((self.numFrames, 3), dtype=float)
        self.kineticEnergy = np.zeros(self.numFrames, dtype=float)
        self.elasticEnergy = np.zeros(self.numFrames, dtype=float)
        self.ljEnergy = np.zeros(self.numFrames, dtype=float)

        self.stage = [{"Name": "", "Frame": [0,0]} for i in range(self.numStages)]
        
        # Initial values
        self.stage[0]["Frame"][0] = 0
        self.stage[-1]["Frame"][1] = self.numFrames

        # Go back and read
        fin.seek(filePos)
        
        # Read every line
        stageIndex = 0
        for i in range(self.numFrames):
            
            # Get line
            line = fin.readline()
            
            # Test for stage
            if "##" in line:
                
                # Get the stage name and limits
                stageName = fin.readline().strip()
                self.stage[stageIndex]["Name"] = stageName
                self.stage[stageIndex]["Frame"][0] = i
                if stageIndex != 0:
                    self.stage[stageIndex-1]["Frame"][1] = i
                stageIndex += 1
                
                # Now get the actual next frame
                fin.readline()
                line = fin.readline()
                
            sline = line.split(",")

            # Assign every value
            hCount = 0
#            Time,CentroidX,CentroidY,CentroidZ,ElasticEnergy,LJEnergy

            for h in headers:
                if h == "Time":
                    self.time[i] = float(sline[hCount])
                elif h == "CentroidX":
                    self.centroid[i][0] = float(sline[hCount])
                elif h == "CentroidY":
                    self.centroid[i][1] = float(sline[hCount])
                elif h == "CentroidZ":
                    self.centroid[i][2] = float(sline[hCount])
                elif h == "LJEnergy":
                    self.ljEnergy[i] = float(sline[hCount])
                elif h == "KineticEnergy":
                    self.kineticEnergy[i] = float(sline[hCount])
                elif h == "ElasticEnergy":
                    self.elasticEnergy[i] = float(sline[hCount])

                hCount += 1
                
        # Close
        fin.close()
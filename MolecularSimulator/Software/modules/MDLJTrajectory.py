#
# These classes read in trajectory data. Each simulation snapshot is a frame, and each frame
# contains position, velocity and force data for each particle. It's a big data structure!
#

import os
import numpy as np

class MDLJTrajectory:
    
    def __init__(self):
        
        self.simulationType = ""
        self.numParticles = 0
        
        self.numFrames = 0
        self.time = None
        self.frame = []

        self.numStages = 0
        self.stage = None

    def load(self, tFname):
        
        # Check filename exists
        if not os.path.exists(tFname):
            print("Cannot read, file does not exist :(")
            raise FileNotFoundError
        
        # Open and check it's the correct file
        fin = open(tFname, "r")

        line = fin.readline().strip()
        if line != "# MS LJGAS TRAJECTORY #":
            print("Cannot read, not the correct filetype :(")
            raise FileNotFoundError

        # Ok it should work from here. Any errors are probably Ben's fault...
        
        # Get simulation type
        self.simulationType = fin.readline().split(",")[1].strip()
            
        # Get num nodes
        self.numParticles = int(fin.readline().split(",")[1].strip())
        
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
            
            elif "#" in line:
                self.numFrames += 1
            
        # Initialise memory structures
        self.time = np.zeros(self.numFrames, dtype=float)
        self.frame = [MDLJTrajectoryFrame(self.numParticles) for i in range(self.numFrames)]

        self.stage = [{"Name": "", "Frame": [0,0]} for i in range(self.numStages)]
        
        # Initial values
        self.stage[0]["Frame"][0] = 0
        self.stage[-1]["Frame"][1] = self.numFrames
        
        # Go back and read
        fin.seek(filePos)

        # Every frame
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

            # This line should be the time
            self.time[i] = float(line)
            
            # Now we need to read the frame information

            # Positions
            line = fin.readline()
            sline = line.split(",")
            count = 0
            for j in range(self.numParticles):
                for k in range(3):
                    self.frame[i].pos[j][k] = float(sline[count])
                    count += 1

            # Velocities
            line = fin.readline()
            sline = line.split(",")
            count = 0
            for j in range(self.numParticles):
                for k in range(3):
                    self.frame[i].vel[j][k] = float(sline[count])
                    count += 1
                    
            # Forces
            line = fin.readline()
            sline = line.split(",")
            count = 0
            for j in range(self.numParticles):
                for k in range(3):
                    self.frame[i].force[j][k] = float(sline[count])
                    count += 1  

            # PBC Wraps
            line = fin.readline()
            sline = line.split(",")
            count = 0
            for j in range(self.numParticles):
                for k in range(3):
                    self.frame[i].wrap[j][k] = int(sline[count])
                    count += 1

            # And the file marker
            fin.readline()

    def unwrap(self):
        
        # For every particle in every frame, we take off all translations caused by the box
        for fIndex in range(self.numFrames):
            for n in range(self.numParticles):
                self.frame[fIndex].pos[n] -= self.frame[fIndex].wrap[n]


class MDLJTrajectoryFrame:
    
    def __init__(self, numParticles):
        
        # Particle data
        self.pos = np.zeros((numParticles,3), dtype=float)
        self.vel = np.zeros((numParticles,3), dtype=float)
        self.force = np.zeros((numParticles,3), dtype=float)
        
        # This is a periodic boundary conditions thing. Learn about it with the SimulationBox
        self.wrap = np.zeros((numParticles,3), dtype=int)


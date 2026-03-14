import numpy as np

#
# This class contains methods to write data to output files
#
class FileWriter:
    
    # The __init__ method runs immediately when a class is instantiated
    # `self.` just means that the variable exists on the object itself,
    # and can be accessed from outside the function with the dot notation
    def __init__(self, tFname, mFname):
        
        self.tFout = None
        self.mFout = None
        
        self.tFname = tFname
        self.mFname = mFname
        
    def openOutputFiles(self, params):
        
        # File pointer objects
        self.tFout = open(self.tFname, "w")
        self.mFout = open(self.mFname, "w")
        
        # Write header
        self.writeTrajectoryHeader(params)
        self.writeMeasurementHeader(params)
        
    def closeOutputFiles(self):
        
        self.tFout.close()
        self.mFout.close()
        
    def writeTrajectoryHeader(self, params):
        
        # Here I am using formatted output to write header data to files
        # '\n' is newline
        self.tFout.write("# MS LJGAS TRAJECTORY #\n")
        
        # '%s' means I'm writing a string. The string is specified in brackets at the end
        # '%d' means I'm writing an integer. The integer is specified in brackets at the end
        self.tFout.write("Simulation Type, %s\n" % (params.simulationType))
        self.tFout.write("numNodes, %d\n" % (params.numParticles))
        
        # Header complete, ready to write simulation frames
        self.tFout.write("###\n")
        
    def writeMeasurementHeader(self, params):
        
        # Here I am using formatted output to write header data to files
        # '\n' is newline
        self.mFout.write("# MS LJGAS MEASUREMENT #\n")

        # Some info about the simulation
        self.mFout.write("Simulation Type, %s\n" % (params.simulationType))

        # Print a titles line which allows us to create a measurement dictionary
        # when reading the data in for analysis
        self.mFout.write("Time,CentroidX,CentroidY,CentroidZ,LJEnergy")
        
        # Only store kinetic energy if inertia exists
        # Saves us storing a load of zeros for no reason
        if params.simulationType == "Langevin":
            self.mFout.write(",KineticEnergy")
        self.mFout.write("\n")
        
        # Header complete, ready to write simulation frames
        self.mFout.write("###\n")
        
    def writeStage(self,header):
        
        # Tags in files where stage is
        self.mFout.write("##\n" + header + "\n##\n")
        self.tFout.write("##\n" + header + "\n##\n")
        
    def writeFrame(self, pSys, box, params, simTime):
        
        # Trajectory first
        self.writeTrajectoryFrame(pSys, params, simTime)
        
        # Then measurement
        self.writeMeasurementFrame(pSys, box, params, simTime)

    def writeTrajectoryFrame(self, pSys, params, simTime):
        
        # Trajectory is just a big list of all positions, velocities and forces
        # Make it so it's easy to read back in
        
        # Time first
        # '%.5f' means I'm writing a float to 5 decimal places 
        self.tFout.write("%.5f\n" % (simTime))
        
        # All positions on a single line: x0,y0,z0,x1,y1,z1...xn,yn,zn
        self.tFout.write("%.5f,%.5f,%.5f" % (pSys.pos[0][0], pSys.pos[0][1], pSys.pos[0][2]))
        for n in range(1, params.numParticles):
            self.tFout.write(",%.5f,%.5f,%.5f" % (pSys.pos[n][0], pSys.pos[n][1], pSys.pos[n][2]))
        self.tFout.write("\n")

        # All velocities on a single line
        self.tFout.write("%.5f,%.5f,%.5f" % (pSys.vel[0][0], pSys.vel[0][1], pSys.vel[0][2]))
        for n in range(1, params.numParticles):
            self.tFout.write(",%.5f,%.5f,%.5f" % (pSys.vel[n][0], pSys.vel[n][1], pSys.vel[n][2]))
        self.tFout.write("\n")
        
        # All forces on a single line
        self.tFout.write("%.5f,%.5f,%.5f" % (pSys.force[0][0], pSys.force[0][1], pSys.force[0][2]))
        for n in range(1, params.numParticles):
            self.tFout.write(",%.5f,%.5f,%.5f" % (pSys.force[n][0], pSys.force[n][1], pSys.force[n][2]))
        self.tFout.write("\n")

        # All wraps on a single line
        self.tFout.write("%d,%d,%d" % (pSys.wrap[0][0], pSys.wrap[0][1], pSys.wrap[0][2]))
        for n in range(1, params.numParticles):
            self.tFout.write(",%d,%d,%d" % (pSys.wrap[n][0], pSys.wrap[n][1], pSys.wrap[n][2]))
        self.tFout.write("\n")
        
        # End of frame marker
        self.tFout.write("#\n")
        
        # Flush forces a write, and clears everything out of the buffers
        self.tFout.flush()
            

    def writeMeasurementFrame(self, pSys, box, params, simTime):
        
        # Measurement should just be 1 line per frame
        centroid = pSys.calculateCentroid()
        ljEnergy = box.getLJEnergy() # This was precalculated in the big box loop to save time

        self.mFout.write("%.5f,%.5f,%.5f,%.5f,%.5f" % (simTime, centroid[0], centroid[1], centroid[2], ljEnergy))
        if params.simulationType == "Langevin":
            ke = pSys.calculateKineticEnergy(params)
            self.mFout.write(",%.2f" % (ke))
        self.mFout.write("\n")
        
        # Flush forces a write, and clears everything out of the buffers
        self.mFout.flush()
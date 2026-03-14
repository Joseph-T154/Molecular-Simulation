#
# This class defines the SimulationRunner
# This purpose of this class is to take in bascially everything as a parameter
# and perform simulation algorithms
#

#
# Module imports
#

# Importing external modules that we'll need
import sys, os
import numpy as np

# Importing the classes we have defined
from modules.SystemParameters import PhysicalConstants
from modules.SystemParameters import SimulationParameters
from modules.PhysicalSystem import PhysicalSystem
from modules.SimulationBox import SimulationBox
from modules.NumericalIntegrator import NumericalIntegrator
from modules.FileIO import FileWriter

class SimulationRunner:
    
    def __init__(self):
        
        # Nothing really needed here, it's just a functional class
        # You could put flags here to see what stage the simulation is at
        # but I'm not going to bother!
        pass

    def runSimulation(self, pSys, params, box, fw, equilibrationTime, numericalIntegrator, thermal=True, lj=True, frameRate = None):

        # Num steps and frame rate corresponding to the equilibration part
        numStepsRequired = int(equilibrationTime / params.dt)

        # Don't print out too many states while equilibrating. Maybe 50 ish?
        if frameRate == None:
            frameRate = int(numStepsRequired / 50.0)
            if frameRate == 0:
                frameRate = 1
        else:
            frameRate = int(frameRate / params.dt)
            
        #
        # Simulation begins here!
        #
        totalSimulationTime = 0.0
        lastWrittenStep = -1
        for step in range(numStepsRequired):

            # Prepare system for a step
            pSys.zeroForces()
            
            if lj:
                box.zeroEnergies()
                
            box.repopulate(pSys)

            #
            # Calculate the net force on the system
            # Thermal needs stochastic (random) force
            # LJ needs Van der Waals LJ interactions
            # Langevin needs drag explicitly adding
            #

            # Thermal noise on each particle
            if thermal:
                pSys.applyThermalForces(params)
                
            # LJ interaction between all pairs of particles
            if lj:
                box.applyLJForces(pSys, params)

            # Drag forces on every particle
            if params.simulationType == "Langevin":
                pSys.applyLocalDragForces(params)

            #
            # Numerically integrate the system
            # 

            # The integrator already knows how to integrate because
            # we told it when we initialised it
            numericalIntegrator.integrate(pSys, params)

            # Do state changes
            
            # Apply boundary conditions
            box.applyBoundaryConditions(pSys, params)

            # Do we output any information? # Writing (to screen or file) is extremely costly
            # so do it relatively sparsely
            if step % frameRate == 0:
                sys.stdout.write("\rStep %d of %d (%.3fns simulated)" % (step + 1, numStepsRequired, totalSimulationTime+params.dt))
                fw.writeFrame(pSys, box, params, totalSimulationTime)
                lastWrittenStep = step

            # Update the simulation time
            totalSimulationTime += params.dt

        # Maybe we need one more write
        if step != lastWrittenStep:
            sys.stdout.write("\rStep %d of %d (%.3fns simulated)" % (step + 1, numStepsRequired, totalSimulationTime+params.dt))
            fw.writeFrame(pSys, box, params, totalSimulationTime)

        # No return value. If I were being a bit more rigid,
        # I would return a pre-defined "Success" variable
        # but that's not good for learning!
        
    
    def runThermalEquilibration(self, pSys, params, box, fw, equilibrationTime, numericalIntegrator):

        # Write equilibration step details to the files
        fw.writeStage("ThermalEquilibration")

        self.runSimulation(pSys, params, box, fw, equilibrationTime, numericalIntegrator, thermal=True, lj=False)
        
    def runLJEquilibration(self, pSys, params, box, fw, equilibrationTime, numericalIntegrator):

        # Write equilibration step details to the files
        fw.writeStage("LJEquilibration")

        self.runSimulation(pSys, params, box, fw, equilibrationTime, numericalIntegrator, thermal=False, lj=True)
        
    def runThermalLJEquilibration(self, pSys, params, box, fw, equilibrationTime, numericalIntegrator):

        # Write equilibration step details to the files
        fw.writeStage("ThermalLJEquilibration")

        self.runSimulation(pSys, params, box, fw, equilibrationTime, numericalIntegrator, thermal=True, lj=True)

    def runProduction(self, pSys, params, box, fw, equilibrationTime, numericalIntegrator, thermal=True, lj=True, frameRate = None):

        # Write equilibration step details to the files
        fw.writeStage("Production")

        self.runSimulation(pSys, params, box, fw, equilibrationTime, numericalIntegrator, thermal=thermal, lj=lj, frameRate = frameRate)
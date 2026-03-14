import os
import numpy as np

#
# This class is simply a container, to hold all the physical constants we'll need
# It might just be Boltzmann's constant, we'll see!
#
class PhysicalConstants:
    
    # The __init__ method runs immediately when a class is instantiated
    # `self.` just means that the variable exists on the object itself,
    # and can be accessed outside the function with the dot notation
    def __init__(self):
        
        self.kb = 1.38e-2

#
# This class is also mostly a container, to hold all the simulation parameters we'll need.
# However, it also has a validate method to check the logical (and physical) consistency
# of the values given
#
class SimulationParameters:
    
    # In this case, the object immediately defines a set of default values for simulation parameters
    # that we will need. In the `validate` function below, these values will be checked for consistency.
    # For some of them, if they do not change from their initial value, the code will register an error.
    def __init__(self):
        
        # Physical data
        self.simulationType = ""
        self.numParticles = 0
        self.volumeFraction = 0.0
        self.temperature = 0.0
        self.viscosity = 0.0
        self.particleDensity = 0.0
        self.particleRadius = 0.0
        
        # Computational data
        self.dt = 0.0
        self.simulationTime = 0.0
        self.frameRate = 0.0
        
        # Box data
        self.boundaryConditions = ""
        self.LJcutoff = 0.0
        self.LJeps = 0.0
        self.LJr0 = 0.0
        
        # Non-primitives
        self.drag = 0.0
        self.kbt = 0.0
        
        self.outputRate = 0
        self.numSteps = 0
        
        # Filenames
        self.tFname_Equilibration = ""
        self.mFname_Equilibration = ""
        self.tFname_Production = ""
        self.mFname_Production = ""
        
    def validate(self):
        
        # Check for logical and physical consistency in parameters
        if self.simulationType != "Brownian" and self.simulationType != "Langevin":
            print("Simulation Type not recognised. We accept 'Brownian' or 'Langevin'")
            return False

        #
        # Filename stuff
        #
        
        # If filenames exist, stop simulation (we could make a backup thing instead in future)
        if os.path.exists(self.mFname_Equilibration):
            print("Equilibration Measurement filename already exists. Please change and try again")
            return False

        # If directory is not writable, stop simulation
        if not os.access(os.path.dirname(os.path.abspath(self.mFname_Equilibration)), os.W_OK):
            print("Equilibration Measurement filename cannot be written to the directory specified. Please change and try again")
            return False

        if os.path.exists(self.tFname_Equilibration):
            print("Equilibration Trajectory filename already exists. Please change and try again")
            return False

        if not os.access(os.path.dirname(os.path.abspath(self.tFname_Equilibration)), os.W_OK):
            print("Equilibration Trajectory filename cannot be written to the directory specified. Please change and try again")
            return False
        
        if os.path.exists(self.mFname_Production):
                print("Production Measurement filename already exists. Please change and try again")
                print("Put this back when you finish debugging!")
#                return False

        if not os.access(os.path.dirname(os.path.abspath(self.mFname_Production)), os.W_OK):
            print("Production Measurement filename cannot be written to the directory specified. Please change and try again")
            return False
        
        if os.path.exists(self.tFname_Production):
            print("Production Trajectory filename already exists. Please change and try again")
            print("Put this back when you finish debugging!")
#            return False

        if not os.access(os.path.dirname(os.path.abspath(self.tFname_Production)), os.W_OK):
            print("Production Trajectory filename cannot be written to the directory specified. Please change and try again")
            return False

        #
        # Physical properties
        #
        if not isinstance(self.numParticles, int):
            print("numParticles must be an integer value")
            return False
            
        if self.numParticles <= 0:
            print("numParticles must be greater than zero")
            return False

        if not isinstance(self.volumeFraction, float) and not isinstance(self.volumeFraction, int):
            print("volumeFraction must be an numerical value")
            return False
            
        if self.volumeFraction <= 0 or self.volumeFraction >= 1:
            print("volumeFraction must be greater than zero and less than 1")
            return False

        if not isinstance(self.temperature, float) and not isinstance(self.temperature, int):
            print("temperature must be an numerical value")
            return False
            
        if self.temperature <= 0:
            print("temperature must be greater than zero")
            return False
            
        if not isinstance(self.viscosity, float) and not isinstance(self.viscosity, int):
            print("viscosity must be an numerical value")
            return False
            
        if self.viscosity <= 0:
            print("viscosity must be greater than zero")
            return False

        if self.particleDensity <= 0:
            print("particleDensity must be greater than zero")
            return False
        
        if not isinstance(self.particleRadius, float) and not isinstance(self.particleRadius, int):
            print("particleRadius must be an numerical value")
            return False
            
        if self.particleRadius <= 0:
            print("particleRadius must be greater than zero")
            return False
            
        # Computational properties
        if not isinstance(self.dt, float) and not isinstance(self.dt, int):
            print("dt must be an numerical value")
            return False
            
        if self.dt <= 0:
            print("dt must be greater than zero")
            return False

        if not isinstance(self.frameRate, float) and not isinstance(self.frameRate, int):
            print("frameRate must be an numerical value")
            return False
            
        if self.frameRate <= 0:
            print("frameRate must be greater than zero")
            return False

        if self.frameRate <= self.dt:
            print("frameRate must be greater than dt")
            return False
            
        if not isinstance(self.simulationTime, float) and not isinstance(self.simulationTime, int):
            print("simulationTime must be an numerical value")
            return False
            
        if self.simulationTime <= 0:
            print("simulationTime must be greater than zero")
            return False

        if self.simulationTime <= self.dt:
            print("simulationTime must be greater than dt")
            return False

        if self.simulationTime <= self.frameRate:
            print("simulationTime must be greater than frameRate")
            return False
         
        #
        # Computational stuff
        #
        
        # Boundary conditions must be a recognised format
        # Change to lower case first
        try:
            self.boundaryConditions = self.boundaryConditions.lower()
        except:
            print("Boundary conditions incorrectly formatted. We need a string!")
            return False
        
        if self.boundaryConditions != "hbc" and self.boundaryConditions != "sbc" and self.boundaryConditions != "pbc":
            print("Boundary conditions must either be hard ('hbc'), soft ('sbc'), or periodic ('pbc')")
            return False
        
        # Box LJ parameters
        if self.LJr0 < 0:
            print("LJr0 must be greater than zero")
            return False
        
        if self.LJeps < 0:
            print("LJeps must be greater than zero")
            return False

        if self.LJcutoff < 0:
            print("LJeps must be greater than zero")
            return False

        if self.LJcutoff <= self.LJr0:
            print("LJeps must be greater than LJr0")
            return False
        
        return True
    
    def setDensity(self, value):
        
        # We are converting to internal units to make the maths work efficiently
        # The conversion is such that acceleration calculated as a=F/m will
        # yield acceleration in nm/ns2 when the force is in pN
        # 1kDa/nm3 = (5/3)e-3 pN.ns2/nm4 (internal units of density)
        # Good practice is to divide away the units coming in, and times them back in going out, so...
#        self.particleDensity = value * (5.0/3.0)*1e-3
        self.particleDensity = value / ((3.0/5.0)*1e3)
                    
    def calculateNonPrimitives(self, pConst):
        
        # Physical
        self.kbt = pConst.kb * self.temperature
        self.drag = 6 * np.pi * self.viscosity * self.particleRadius  # Stokes drag: spherical particles, no hydrodynamic interactions
        self.particleMass = self.particleDensity * (4.0/3.0) * np.pi * self.particleRadius**3
        
        # Computational
        self.numSteps = int(np.ceil(self.simulationTime / self.dt)) # Total number of simulation iterations required
        self.outputRate = int(np.ceil(self.frameRate / self.dt)) # Output every 'outputRate' simulation steps
        
        self.simulationTime = self.dt * self.numSteps  # Recalculate to account for the int casting
        self.frameRate = self.dt * self.outputRate

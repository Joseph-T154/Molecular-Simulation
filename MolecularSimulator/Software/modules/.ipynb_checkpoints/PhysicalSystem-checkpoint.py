import numpy as np
from modules.PhysicalObjects import LinearSpring

class PhysicalSystem:

    # The __init__ method runs immediately when a class is instantiated
    # `self.` just means that the variable exists on the object itself,
    # and can be accessed outside the function with the dot notation
    def __init__(self, params):
        
        # Store numParticles on object for fast lookup
        self.numParticles = params.numParticles
        
        # Particle data
        self.pos = np.zeros((self.numParticles,3), dtype=float)
        self.vel = np.zeros((self.numParticles,3), dtype=float)
        self.force = np.zeros((self.numParticles,3), dtype=float)

        # This is a periodic boundary conditions thing. Learn about it with the SimulationBox
        self.wrap = np.zeros((self.numParticles,3), dtype=int)

        # Store numSprings (and numPolymers) for fast lookup
        self.numSprings = params.numSprings
        self.numPolymers = self.numParticles - self.numSprings
        self.numMonomers = int(self.numParticles / self.numPolymers)
        
        # Spring data (physical parameters and connectivity)
        self.spring = [LinearSpring(params.polymerK, params.polymerL) for i in range(self.numSprings)]

        sIndex = 0
        for i in range(self.numPolymers):
            for j in range(self.numMonomers -  1):
                a = self.numMonomers * i + j
                b = self.numMonomers * i + j + 1
                self.spring[sIndex].setIndices(a,b)
                sIndex += 1
                
        self.elasticEnergy = 0.0
        
    def zeroForces(self):
        self.force.fill(0.0)

    def zeroEnergies(self):
        self.elasticEnergy = 0.0
        for s in range(self.numSprings):
            self.spring[s].elasticEnergy = 0.0
            
    # This method applies thermal forces to each particle
    # It utilises the fluctuation-dissipation theorem
    def applyThermalForces(self, params):
    
        # All particles are independent in these simulations. So, loop over particles in isolation 
        # and apply stochastic force component in each direction
        for n in range(self.numParticles):

            for j in range(3):

                # Thermal force (as per fluctuation-dissipation theorem) is a gaussian random
                # with mean 0, variance 2*kbt*drag/dt
                # This results in 0.5kT of energy (on average) per degree of freedom
                self.force[n][j] += np.random.normal(0, np.sqrt(2*params.kbt*params.drag / params.dt))

    #
    # This method applies elastic forces to all particles connected to elastic objects
    #
    def applyElasticForces(self, params):
        
        for s in range(self.numSprings):

            # Get spring state
            st = self.spring[s].currentState
            
            # Get extension (and consider the PBC wrap)
            n = self.spring[s].n
            rij = (self.pos[n[1]] - self.wrap[n[1]]) - (self.pos[n[0]] - self.wrap[n[0]])
            r = np.linalg.norm(rij)
            rijHat = rij / r
            
            # Force magnitude
            dr = r - (self.spring[s].l[st] + 2*params.particleRadius)
            fMag = self.spring[s].k[st] * dr
            
            # Apply forces in correct direction
            force = fMag * rijHat
            self.force[n[0]] += force
            self.force[n[1]] -= force

            # Get energy for measurement files
            self.spring[s].elasticEnergy = 0.5 * self.spring[s].k[st] * dr**2
            self.elasticEnergy += self.spring[s].elasticEnergy

        # Finalise elastic energy
        self.elasticEnergy *= 0.5

    # This method applies local drag forces to each particle
    # This function is only needed if we are using Langevin (inertial) integration
    # This function is not needed if we are performing generalised hydrodynamics
    def applyLocalDragForces(self, params):
    
        # All particles are independent in these simulations. So, loop over particles in isolation 
        # and apply velocity-dependent drag force
        for n in range(self.numParticles):

            self.force[n] -= params.drag * self.vel[n]
            
    def calculateCentroid(self):
        
        # This is centre of position
        # It is NOT centre of mass
        centroid = np.zeros(3,dtype=float)
        for n in range(self.numParticles):

            centroid += self.pos[n]
            
        return centroid / self.numParticles
    
    def calculateKineticEnergy(self, params):
        
        ke = 0.0
        for n in range(self.numParticles):

            # If we had a simulation where the particles hadd different masses,
            # the mass value would have to come inside the summation
            ke += np.dot(self.vel[n], self.vel[n])
            
        ke *= 0.5 * params.particleMass
        return ke

    def getElasticEnergy(self):

        return self.elasticEnergy

    def stateChange(self, params):

        # For all springs in the system, get the energy and maybe change state
        for s in range(self.numSprings):

            Ei = self.spring[s].elasticEnergy
            i = self.spring[s].currentState
            n = self.spring[s].n
            
            rij = (self.pos[n[1]] - self.wrap[n[1]]) - (self.pos[n[0]] - self.wrap[n[0]])
            r = np.linalg.norm(rij)
            dr = r - (self.spring[s].l[i] + 2*params.particleRadius)

            tProbsNow = np.zeros(params.numStates, dtype=float)
            
            # Get all transition probabilities based on the current energy
            totalProb = 0.0
            for j in range(params.numStates):

                if i == j:
                    continue

                Ej = 0.5 * self.spring[s].k[j] * dr**2
                tProbsNow[j] = params.pureTProb[i][j] * np.exp(params.fij*(Ei-Ej) / params.kbt)
                totalProb += tProbsNow[j]

            # Time resolution
#            print(tProbsNow)
            if totalProb > 1.0:
                print("Sum of instantaneous transition probabilities greater than 1. Lower your timestep or lower your rates to resolve these interactions.")
                raise RuntimeError

            # Finalise probabilities
            tProbsNow[i] = 1 - totalProb

            # Change based on random number
            rTest = np.random.rand()

            totalProb = 0.0
            for j in range(params.numStates):
                totalProb += tProbsNow[j]
                if rTest <= totalProb:
                    self.spring[s].currentState = j
                    break
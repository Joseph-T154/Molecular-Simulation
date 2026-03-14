import numpy as np

class PhysicalSystem:

    # The __init__ method runs immediately when a class is instantiated
    # `self.` just means that the variable exists on the object itself,
    # and can be accessed outside the function with the dot notation
    def __init__(self, numParticles):
        
        # Store numParticles on object for fast lookup
        self.numParticles = numParticles
        
        # Particle data
        self.pos = np.zeros((numParticles,3), dtype=float)
        self.vel = np.zeros((numParticles,3), dtype=float)
        self.force = np.zeros((numParticles,3), dtype=float)

        
        # This is a periodic boundary conditions thing. Learn about it with the SimulationBox
        self.wrap = np.zeros((numParticles,3), dtype=int)

    def zeroForces(self):
        self.force.fill(0.0)
    
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
#
# This class exists to perform integration operators
# Numerical integration is quite abstract, and so can be defined
# separately from the PhysicalSystem itself
#

class NumericalIntegrator:
    
    def __init__(self, integrationType):
        
        self.integrate = None
        
        # Bind the correct function
        if integrationType == "Brownian":
            self.integrate = self.brownianIntegrate
        elif integrationType == "Langevin":
            self.integrate = self.langevinIntegrate
        else:
            
            # This condition should've been caught before it reaches here, but you never know
            raise TypeError("Unrecognised simulation type. Cannot integrate")
            
    def brownianIntegrate(self, physicalSystem, parameters):

        # Integrate with respect to drag. No inertia!
        for n in range(physicalSystem.numParticles):

            # With no inertia, the velocity is not persistent, is set directly from the force
            # Hence, kinetic energy is not properly defined

            # If particles had different radii, drag would be per particle, not global
            # v(t+dt) = f(t) / drag
            physicalSystem.vel[n] = physicalSystem.force[n] / parameters.drag

            # Position is dynamically updated rather than set. This is a core difference
            # between molecular dynamics and Monte Carlo simulations
            # x(t+dt) = x(t) + v(t+dt)dt
            physicalSystem.pos[n] += physicalSystem.vel[n] * parameters.dt

    def langevinIntegrate(self, physicalSystem, parameters):

        # Integrate with respect to mass
        for n in range(physicalSystem.numParticles):

            # With inertia, the velocity is persistent, and it updated from the force
            # Hence, kinetic energy is properly defined

            # If particles had different radii, drag would be per particle, not global
            # v(t+dt) = v(t) * (f(t) / m)*dt
            acc = physicalSystem.force[n] / parameters.particleMass
            physicalSystem.vel[n] += acc * parameters.dt

            # Position is dynamically updated rather than set. This is a core difference
            # between molecular dynamics and Monte Carlo simulations
            # x(t+dt) = x(t) + v(t+dt)dt
            physicalSystem.pos[n] += physicalSystem.vel[n] * parameters.dt
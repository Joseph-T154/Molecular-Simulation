#
# This class defines the SimulationBox
# The function of a simulation box is to optimise physics calculations
# and create boundary conditions. Boundary conditions have physical effects
# which must be taken into account! Optimisations, on the other hand, do not.
# All they do is speed up the calculations, they do not change them in any way.
# For example, they can be used to ignore negligible forces.
#

import numpy as np
from enum import Enum

# This class is an enumerated collection, inherited from the Enum class
# I am basically setting variables so I can use them later
# Ask if you're interested and don't know what class inheritance is!
class BoundaryConditions(Enum):
    
    NONE = 0
    HBC = 1
    SBC = 2
    PBC = 3
    
class Voxel:

    # The __init__ method runs immediately when a class is instantiated
    # `self.` just means that the variable exists on the object itself,
    # and can be accessed outside the function with the dot notation
    def __init__(self):
        
        # A voxel should know it's own structure indices for fast lookup
        self.indices = np.zeros(3, dtype=int)
        self.index = 0
        
        # A voxel needs to contain something
        self.numObjects = 0
        self.nodeIndex = set()
        
        # A voxel needs to know what it's connected to
        # This will change depending on boundary conditions
        self.numNeighbours = 0
        self.neighbours = set()
        
        # If we loop over all voxels, then all neighbours, we will double count
        # To address this, we build a balanced set of unique neighbours for pair-pair calculations
        # Thus, if voxel i does calculations with j, j will not do them with i
        self.numCalculationNeighbours = 0
        self.calculationNeighbours = set()
    
        # Is this voxel an edge voxel?
        self.edge = False
        
    def clear(self):
        
        # Sets have a clear function! Nice
        self.nodeIndex.clear()
        
    def connectToVoxel(self, index):
        self.neighbours.add(index)
        self.numNeighbours += 1

    def connectToCalculationVoxel(self, index):
        self.calculationNeighbours.add(index)
        self.numCalculationNeighbours += 1

    def isConnected(self, index):
        return index in self.neighbours
    
    def isCalculationConnected(self, index):
        return index in self.calculationNeighbours
    
    def addObject(self, index):
        
        # Add particle to the list
        self.nodeIndex.add(index)
        self.numObjects += 1
    
class SimulationBox:
    
    # The __init__ method runs immediately when a class is instantiated
    # `self.` just means that the variable exists on the object itself,
    # and can be accessed outside the function with the dot notation
    def __init__(self, bc, length=0.0, resolution=0.0):
        
        # Some box quantities
        self.LJenergy = 0.0
        
        # Get the boundary conditions
        if bc == "hbc":
            self.bc = BoundaryConditions.HBC
        elif bc == "sbc":
            self.bc = BoundaryConditions.SBC
        elif bc == "pbc":
            self.bc = BoundaryConditions.PBC
        else:
            raise TypeError("Unrecognised boundary conditions: %s" % (bc))
            
        # Check input parameters
        eps = 0.01  # This is a lower size limit, chosen arbitrarily
        if length < eps or resolution < eps:
            raise ValueError("Box length or resolution is too small. Please try again")
            
        if length < resolution:
            raise ValueError("Box length must be greater than the resolution. Please try again")
        
        # Calculate the number of voxels (3D pixels) required in the simulation box
        numVoxels1D = int(np.ceil(length/resolution))   # Get the integer limit of the number of voxels
        
        # Check this, as it should be at least 3 or greater for pbc (so a voxel is not connected to itself)
        if self.bc == BoundaryConditions.PBC and numVoxels1D < 3:
            raise ValueError("Insufficient number of voxels for PBC. Make resolution smaller.")
            
        self.dimension = np.array([length,length,length], dtype=float)
        self.volume = self.dimension[0] * self.dimension[1] * self.dimension[2]
        self.numVoxels = np.array([numVoxels1D,numVoxels1D,numVoxels1D], dtype=int)
        self.totalNumVoxels =  int(numVoxels1D**3)

        self.resolution = length / numVoxels1D  # Recalculate to account for integer rounding above

        #
        # Build the voxels
        #
        
        # Voxels are basically sets that will contain particle indices
        # They represent 3D boxes containing particles
        # The populate method will take in the pSys and refill this object
        self.buildVoxels()

        # Now we need to build a lookup table for when we do calculations, to speed them up
        if self.bc == BoundaryConditions.PBC:
            self.buildImageLookup()
        
    def buildVoxels(self):

        # Create the objects
        self.voxel = [Voxel() for i in range(self.totalNumVoxels)]
        self.imageTranslation = {}
        
        # Connect them together based on the boundary conditions
        edge = False
        for i in range(self.numVoxels[0]):
            
            # Is it an edge?
            if i == 0 or i == self.numVoxels[0] - 1:
                edge = True
                
            for j in range(self.numVoxels[1]):
                
                # Is it an edge?
                if j == 0 or j == self.numVoxels[1] - 1:
                    edge = True
                    
                for k in range(self.numVoxels[2]):

                    # Is it an edge?
                    if k == 0 or k == self.numVoxels[2] - 1:
                        edge = True
                        
                    # Get the list index
                    fullIndex = self.indexTransform3To1(i,j,k)
                    
                    # Set the voxel indices
                    self.voxel[fullIndex].index = fullIndex
                    self.voxel[fullIndex].indices[0] = i
                    self.voxel[fullIndex].indices[1] = j
                    self.voxel[fullIndex].indices[2] = k
                    
                    # And its edge flag
                    self.voxel[fullIndex].edge = edge
                    
                    # What is this voxel connected to? Everything one away in every direction!
                    # Incidentally, you should never have 6 embedded loops except in extreme conditions
                    # and where you understand the limits of each loop!
                    # In this case, the "buildBox" operation happens only once in the simulation
                    for dx in range(-1,2,1):
                        
                        # Other voxel indices
                        i2 = i + dx
                        
                        # In each loop, check boundary conditions first!     
                        # Boundary conditions. What do we do if index goes past end of box?
                        # If periodic, loop around box. If not, don't connect
                        if self.bc == BoundaryConditions.HBC:
                            if i2 >= self.numVoxels[0] or i2 < 0:
                                continue
                        elif self.bc == BoundaryConditions.PBC:
                            if i2 >= self.numVoxels[0]:
                                i2 = (i2 + self.numVoxels[0]) % self.numVoxels[0]

                            elif i2 < 0:
                                i2 = (i2 + self.numVoxels[0]) % self.numVoxels[0]

                        for dy in range(-1,2,1):

                            # Other voxel indices
                            j2 = j + dy

                            # Boundary conditions
                            if self.bc == BoundaryConditions.HBC:
                                if j2 >= self.numVoxels[1] or j2 < 0:
                                    continue
                            elif self.bc == BoundaryConditions.PBC:
                                if j2 >= self.numVoxels[1]:
                                    j2 = (j2 + self.numVoxels[1]) % self.numVoxels[1]

                                elif j2 < 0:
                                    j2 = (j2 + self.numVoxels[1]) % self.numVoxels[1]

                            for dz in range(-1,2,1):
                                
                                # Other voxel indices
                                k2 = k + dz

                                # Boundary conditions
                                if self.bc == BoundaryConditions.HBC:
                                    if k2 >= self.numVoxels[2] or k2 < 0:
                                        continue
                                elif self.bc == BoundaryConditions.PBC:
                                    if k2 >= self.numVoxels[2]:
                                        k2 = (k2 + self.numVoxels[2]) % self.numVoxels[2]
                                    elif k2 < 0:
                                        k2 = (k2 + self.numVoxels[2]) % self.numVoxels[2]
                                        
                                # Finally, do not connect the voxel to itself
                                if dx == 0 and dy == 0 and dz == 0:
                                    continue
                                    
                                # If we get here, we can connect!
                                fullNeighbourIndex = self.indexTransform3To1(i2,j2,k2)

                                # Add both to each others neighbour list
                                # Voxel connection objects are sets, so no need to worry
                                # about double-counting
                                self.voxel[fullIndex].connectToVoxel(fullNeighbourIndex)
                                self.voxel[fullNeighbourIndex].connectToVoxel(fullIndex)
                                
                                # Now, sort out calculationNeighbours. One voxel, and only one voxel,
                                # should connect to the other
                                
                                # Is it already connected though? No double counting here!
                                if not self.voxel[fullIndex].isCalculationConnected(fullNeighbourIndex) and not self.voxel[fullNeighbourIndex].isCalculationConnected(fullIndex):

                                    # Randomise which way it is to keep it roughly balanced
                                    testVal = np.random.randint(0,2)
                                    if testVal == 0:
                                        self.voxel[fullIndex].connectToCalculationVoxel(fullNeighbourIndex)
                                    else:
                                        self.voxel[fullNeighbourIndex].connectToCalculationVoxel(fullIndex)

    def buildImageLookup(self):
        
        # Get the lookup
        self.imageTranslation = {}
        
        # Loop over all voxels and their neighbours
        for cV in self.voxel:

            # Now do all interactions with the neighbours
            for nIndex in cV.calculationNeighbours:

                # Get the neighbouring voxel
                nV = self.voxel[nIndex]
                
                # If they are both edges, we need to add to the lookup table
                if cV.edge and nV.edge:
                    
                    # Initialise
                    translation = np.zeros(3,dtype=float)
                    
                    # Check each direction
                    for i in range(3):
                        if cV.indices[i] == 0 and nV.indices[i] == self.numVoxels[i] - 1:
                            translation[i] -= self.dimension[i]
                        elif cV.indices[i] == self.numVoxels[i] - 1 and nV.indices[i] == 0:
                            translation[i] += self.dimension[i]

                    # Add
                    self.imageTranslation[(cV.index,nV.index)] = translation

                        
                
    def clear(self):
        
        # Empty all voxel contents. Do not change connectivity!
        for i in range(self.totalNumVoxels):
            self.voxel[i].clear()
        
    def repopulate(self, pSys):
        
        # Clear, then populate
        self.clear()
        self.populate(pSys)
        
    def populate(self, pSys):
        
        # This method will assign all particle indices to the box
        # It will not do anything about boundary conditions, that is separate
        # Hence, if a particle is outside the box...it shouldn't be!
        # We raise an error
        
        # Reset the box
        self.clear()
        
        # Loop over all particles in the system
        for n in range(pSys.numParticles):
           
            # Which voxel should it be in, in x,y,z index co-ordinates
            
            # Get index via division operation
            indices = [int(vIndex) for vIndex in pSys.pos[n] / self.resolution]
            fullIndex = self.indexTransform3To1(indices[0], indices[1], indices[2])
            
            # Add the index to the voxel
            try:
                self.voxel[fullIndex].addObject(n)
            except:
                print(pSys.pos[n], self.resolution, indices, fullIndex,self.dimension)
                raise
                
    def indexTransform3To1(self, x,y,z):
        
        # The box is such that we count all z's, then all y's, then all x's
        return (x * self.numVoxels[1] + y) * self.numVoxels[2] + z
    
    def indexTransform1To3(self, index):

        # Consider the above method. How do we reverse it?
        indices = [None, None, None]
        
        # z is the remainder
        indices[2] = index % self.numVoxels[2]
        
        # y is the remainder once all the z has been taken off
        indices[1] = ((index - indices[2]) / self.numVoxels[2]) % self.numVoxels[1];
        
        # x is the remainder once all the y and z have been taken off        
        indices[0] = (((index - indices[2]) / self.numVoxels[2]) - indices[1]) / self.numVoxels[1];

        return indices;

    def randomInitialise(self, pSys, params):
        
        # This method will randomly distribute the particles throughout the box
        # It will also make sure none of the radii cross the box boundary
        # in case of non-periodic boundary conditions
        
        # For every node, put in box somewhere
        for n in range(pSys.numParticles):
            
            # In every direction
            for i in range(3):

                # Random number between radius and boxLength - radius
                # Random * scale + offset
                pSys.pos[n][i] = np.random.rand() * (self.dimension[i] - 2.0 * params.particleRadius) + params.particleRadius
                
    def centerInitialise(self, pSys, params):
        
        # This method will put every single particle right at the centre of the box
        # Extremely high entropy state
        for n in range(pSys.numParticles):
            pSys.pos[n] = self.dimension / 2.0
                
    def uniformInitialise(self, pSys, params):
        
        # This method will uniformly distribute the particles throughout the box
        # It will also make sure none of the radii cross the box boundary
        # in case of non-periodic boundary conditions
        
        # Cube root the number of particles and round up
        numX = int(np.ceil(pSys.numParticles**(1.0/3.0)))
        steps = (self.dimension - 2.0 * params.particleRadius) / numX
        aPos = np.zeros(3, dtype=float)
        
        # For every node, put in box somewhere
        numPlaced = 0
        
        # Reset the i run
        aPos[0] = params.particleRadius
        for i in range(numX):
            
            # Are we done?
            if numPlaced == pSys.numParticles:
                break
                
            # Reset the j run
            aPos[1] = params.particleRadius
            
            for j in range(numX):
                
                # Are we done?
                if numPlaced == pSys.numParticles:
                    break

                # Reset the k run
                aPos[2] = params.particleRadius
                for k in range(numX):
                    
                    # Are we done?
                    if numPlaced == pSys.numParticles:
                        break
                        
                    # Place particle
                    pSys.pos[numPlaced] = aPos
                    
                    # Increment
                    aPos[2] += steps[2]
                    numPlaced += 1

                # Increment
                aPos[1] += steps[1]
        
            # Increment
            aPos[0] += steps[0]

    # This method could use function binding when boundary conditions are set, but I want to be 
    # clear exactly what is happening
    def applyBoundaryConditions(self, pSys, params):
        
        # Different method depending on the boundary conditions
        if self.bc == BoundaryConditions.HBC:
            self.applyHardBoundaryConditions(pSys, params)
        elif self.bc == BoundaryConditions.SBC:
            self.applySoftBoundaryConditions(pSys, params)
        elif self.bc == BoundaryConditions.PBC:
            self.applyPeriodicBoundaryConditions(pSys, params)
            
    def applyHardBoundaryConditions(self, pSys, params):
        
        # Hard walls conserve momentum perfectly and reflect the velocity
        # For perfect accuracy, we would reverse the velocity path and work
        # out exactly where the particle would have gone. For speed, with small dt,
        # we will just move the particle back in the box
        
        # Ideally, we would only need to check edge voxels, but just in case, we will
        # check all particles
        for n in range(pSys.numParticles):
            
            # And directions
            for i in range(3):
                
                if pSys.pos[n][i] + params.particleRadius > self.dimension[i]:
                    pSys.pos[n][i] = self.dimension[i] -  params.particleRadius
                    pSys.vel[n][i] *= -1
                elif pSys.pos[n][i] -  params.particleRadius < 0:
                    pSys.pos[n][i] = params.particleRadius
                    pSys.vel[n][i] *= -1                    

    def applySoftBoundaryConditions(self, pSys, params):
        
        # Soft walls absorb momentum perfectly and set velocity to zero
        # For perfect accuracy, we would reverse the velocity path and work
        # out exactly where the particle would have stopped. For speed, with small dt,
        # we will just move the particle back in the box
        
        # Ideally, we would only need to check edge voxels, but just in case, we will
        # check all particles
        for n in range(pSys.numParticles):
            
            # And directions
            for i in range(3):
                
                if pSys.pos[n][i] + params.particleRadius > self.dimension[i]:
                    pSys.pos[n][i] = self.dimension[i] - params.particleRadius
                    pSys.vel[n][i] = 0.0
                elif pSys.pos[n][i] - params.particleRadius < 0:
                    pSys.pos[n][i] = params.particleRadius
                    pSys.vel[n][i] = 0.0
                    
    def applyPeriodicBoundaryConditions(self, pSys, params):
        
        # Periodic walls tranport the particlearound the box. Momentum conserved by default
        
        # Ideally, we would only need to check edge voxels, but just in case, we will
        # check all particles
        for n in range(pSys.numParticles):
            
            # And directions
            for i in range(3):
                
                if pSys.pos[n][i] >= self.dimension[i]:
                    pSys.pos[n][i] -= self.dimension[i]
                    
                    # Store the number of times the box has gone round
                    pSys.wrap[n][i] -= self.dimension[i]
                    
                elif pSys.pos[n][i] < 0:
                    pSys.pos[n][i] += self.dimension[i]

                    # Store the number of times the box has gone round
                    pSys.wrap[n][i] += self.dimension[i]
                    
    
    def applyLJForces(self, pSys, params):
        
        #
        # This is the method that will take the most time in the entire simulation
        # I guarantee it. Therefore, we will do whatever is needed to speed it up
        #
        
        # Setup a loop over all unique pairs (no double counting)
        # This should be a loop over all voxels, and a subsequent loop over
        # unique neighbours
        for cV in self.voxel:
        
            # Do all interactions within the cV
            for cpIndex1 in cV.nodeIndex:
                    for cpIndex2 in cV.nodeIndex:
                        
                        # Skip if same
                        if cpIndex1 == cpIndex2:
                            continue
                        
                        # LJ if not same
                        self.doLJInteraction(cpIndex1, cpIndex2, pSys, params, translation = (None,None))

            # Now do all interactions with the neighbours
            for nIndex in cV.calculationNeighbours:

                # Get the neighbouring voxel
                nV = self.voxel[nIndex]

                # Now, double loop over all particles

                # If both of these are edge voxels, and it's PBC, we have to work out the translation image
                if self.bc == BoundaryConditions.PBC:

                    # If they are both edge voxels, get the translation
                    if cV.edge and nV.edge:
                        translation = self.imageTranslation[(cV.index,nV.index)]
                        
                    for cpIndex1 in cV.nodeIndex:
                            for npIndex2 in nV.nodeIndex:
                                self.doLJInteraction(cpIndex1, npIndex2, pSys, params, translation = (1,translation))
                else:
                    for cpIndex1 in cV.nodeIndex:
                            for npIndex2 in nV.nodeIndex:
                                self.doLJInteraction(cpIndex1, npIndex2, pSys, params, translation = (None,None))
   
                                
    def doLJInteraction(self, cpIndex1, cpIndex2, pSys, params, translation = (None,None)):
        
        # Two particles, image translation for PBC
        # Distance is the surface to surface distance
        # I have removed the hard repulsive bit and replaced with a harmonic potential for simulation stability
        
        # Calculate the vector separation
        if translation[0] == None:
            rcc = pSys.pos[cpIndex2] - pSys.pos[cpIndex1]
        else:
            rcc = pSys.pos[cpIndex2] - pSys.pos[cpIndex1] + translation[1]
        
        # Magnitudes
        rMagcc = np.linalg.norm(rcc)
        rMagss = rMagcc - 2 * params.particleRadius
        
        # Check against cutoff distance (optimisation to save more time)
        if rMagss > params.LJcutoff:
            return
        
        # Unit vector
        try:
            rHat = rcc / rMagcc
        except ZeroDivisionError:
            
            # Random direction then
            rHat = np.random.rand(3)
            rHat /= np.linalg.norm(rHat)
            
        # Full lennard-jones or harmonic approximation?
        if rMagss > params.LJr0:
            
            # Ok, now lets calculate it!
            r0Onr = params.LJr0 / rMagss
            r0Onr6 = r0Onr**6

            # Energy first. Add to the total (as we will be summing over all interactions)
            self.LJenergy += params.LJeps * r0Onr6 * (r0Onr6 - 2)

            # Forces
            forceMag = 12 * (params.LJeps / rMagss) * r0Onr6 * (r0Onr6 - 1)
            pSys.force[cpIndex2] += forceMag * rHat
            pSys.force[cpIndex1] -= forceMag * rHat

        else:
            
            # Harmonic approximation
            LJSpringConstant = 72 * params.LJeps / (params.LJr0**2)
            self.LJenergy += 0.5 * LJSpringConstant * (rMagss - params.LJr0)**2 - params.LJeps
            forceMag = -LJSpringConstant * (rMagss - params.LJr0)
            pSys.force[cpIndex2] += forceMag * rHat
            pSys.force[cpIndex1] -= forceMag * rHat
        
    def getLJEnergy(self):
        
        return self.LJenergy
    
    def zeroEnergies(self):
        
        self.LJenergy = 0.0
#
# This class defines the an object which can test the simulation at various stages as it runs
# It will not give us a "Yes" or "No" as to whether we are good to continue
# Rather, it will give us metrics and graphs which we can use to determine whether to continue
# Due to the computational nature of simulation, and the inherent potential for numerical inaccuracy
# simulation is as much an art as it is a science. This isn't a bad thing. If anything, all of science
# is a much art as it is science! We do our best to explain things as objectively as we can :)
#

import sys, os
import numpy as np
from matplotlib import pyplot as plt

from modules.MDLJMeasurement import MDLJMeasurement
from modules.MDLJTrajectory import MDLJTrajectory

# These are global parameters for matploblib, so I'll put them outside the class
SSMALL = 6
SMALL = 8
MEDIUM = 10

def maxwellBoltzmann(v, params):
    
    preFactor = (params.particleMass / (2.0*np.pi*params.kbt))**(3.0/2.0)
    preFactorSolidAngle = preFactor * 4.0 * np.pi * v**2
    prob = preFactorSolidAngle * np.exp(-(params.particleMass / (2.0*params.kbt)) * v**2)
    
    # Return frequency, not probability
    return prob# * params.numParticles
    
class SimulationTester():
    
    def __init__(self):
        pass
    
    def testThermalEquilibration(self, pSys, params):

        # We will test the equilibration by measuring the spatial distribution of particles
        # Should be uniform!
        sDistFig, sDistaxs = plt.subplots(3, 1, figsize=(3.2,3.2))
        vDistFig = plt.figure(figsize=(3.2,3.2))
        kineticEnergyFig = plt.figure(figsize=(3.2,3.2))
        
        # Distribution figure first
        plt.figure(sDistFig.number)
        
        # Get transposed array data
        tPos = np.transpose(pSys.pos)

        # Plot for each direction
        directionLabel = ["x", "y", "z"]
        maxFreq = 0
        for i in range(3):

            # Plot histogram
            axisData = sDistaxs[i].hist(tPos[i])
            numBins = axisData[1].size - 1

            # Plot expected (uniform)
            xVals = np.linspace(axisData[1][0],axisData[1][-1],100)
            yVals = np.ones(xVals.size) * params.numParticles / numBins
            sDistaxs[i].plot(xVals, yVals, "--k", label="Expected")

            # Prettify
            sDistaxs[i].set_xlabel("%s-Position (nm)" % (directionLabel[i]), fontsize=SMALL)
            sDistaxs[i].set_ylabel("Freq", fontsize=SMALL)
            sDistaxs[i].legend(loc=0, fontsize=SSMALL)

            # Get y-axis for maxFreq (uniform plots)
            yMax = sDistaxs[i].get_ylim()[1]
            if yMax > maxFreq:
                maxFreq = yMax
        
        for i in range(3):
            sDistaxs[i].set_ylim(0,maxFreq)
            plt.tight_layout()

        oFname = os.path.splitext(params.mFname_Equilibration)[0] + "_ThermalEquilPosDist.png"
        plt.savefig(oFname, dpi=100)
        
        # Also, for Langevin, we can measure speed distribution
        # Should be Maxwell-Boltzmann
        if params.simulationType == "Langevin":

            # Vel figure
            plt.figure(vDistFig.number)

            # Get speeds
            speedArray = np.array([np.linalg.norm(v) for v in pSys.vel])

            # Plot distribution
            plt.hist(speedArray, density=True)
            
            # Plot expected (maxwell-boltzmann)
            xVals = np.linspace(0,plt.xlim()[1],100)
            yVals = maxwellBoltzmann(xVals, params)
            plt.plot(xVals, yVals, "--k", label="Maxwell-Boltzmann")
            
            # Prettify
            plt.xlabel("Speed (nm / ns)", fontsize=MEDIUM)
            plt.ylabel("Probability", fontsize=MEDIUM)
            plt.legend(loc=0, fontsize=SMALL)
            
            oFname = os.path.splitext(params.mFname_Equilibration)[0] + "_ThermalEquilVelDist.png"
            plt.savefig(oFname, dpi=100)
            
            # Kinetic energy figure
            plt.figure(kineticEnergyFig.number)
            
            # We want kinetic energy over time, so for this one, we'll need to open the external file
            meas = MDLJMeasurement()
            meas.load(params.mFname_Equilibration)
            
            # Get the Thermal equilibration limits
            start = -1
            end = -1
            for stage in meas.stage:
                if stage["Name"] == "ThermalEquilibration":
                    start = stage["Frame"][0]
                    end = stage["Frame"][1]
            
            if start != -1:
                
                # Plot measured
                plt.plot(meas.time[start:end], meas.kineticEnergy[start:end], "-bx", label="Kinetic Energy")

                # Plot expected (half kT per degree of freedom)
                xVals = np.linspace(0,plt.xlim()[1],100)
                yVals = 1.5 * params.kbt * params.numParticles * np.ones(100, dtype="float")
                plt.plot(xVals, yVals, "--k", label="Equipartition")

                # Prettify
                plt.xlabel("Time (ns)", fontsize=MEDIUM)
                plt.ylabel("Kinetic Energy (pN.nm)", fontsize=MEDIUM)
                plt.legend(loc=0, fontsize=SMALL)

                oFname = os.path.splitext(params.mFname_Equilibration)[0] + "_ThermalEquilKineticEnergy.png"
                plt.savefig(oFname, dpi=100)            
            
        plt.show()
        plt.clf()
        
    def testLJEquilibration(self, pSys, params):

        # For this one, we just need the LJ energy profile
        ljEnergyFig = plt.figure(figsize=(3.2,3.2))

        # We want lj energy over time, so for this one, we'll need to open the external file
        meas = MDLJMeasurement()
        meas.load(params.mFname_Equilibration)

        # Get the LJ equilibration limits
        start = -1
        end = -1
        for stage in meas.stage:
            if stage["Name"] == "LJEquilibration":
                start = stage["Frame"][0]
                end = stage["Frame"][1]

        if start != -1:

            # Plot measured
            plt.plot(meas.time[start:end], meas.ljEnergy[start:end], "-bx", label="LJ Energy")

            # Expected unknown, but should be stable
 
            # Prettify
            plt.xlabel("Time (ns)", fontsize=MEDIUM)
            plt.ylabel("LJ Energy (pN.nm)", fontsize=MEDIUM)
            plt.legend(loc=0, fontsize=SMALL)

            oFname = os.path.splitext(params.mFname_Equilibration)[0] + "_LJEquilLJEnergy.png"
            plt.savefig(oFname, dpi=100)           

        plt.show()
        plt.clf()

    def testThermalLJEquilibration(self, pSys, params):

        # For this one, test all energies and particle distributions
        kineticEnergyFig = plt.figure(figsize=(3.2,3.2))
        ljEnergyFig = plt.figure(figsize=(3.2,3.2))
 
        # We want lj energy over time at least, maybe kinetic too so for this one, we'll need to open the external file
        meas = MDLJMeasurement()
        meas.load(params.mFname_Equilibration)
                
        # Also, for Langevin, we can measure speed distribution
        # Should be Maxwell-Boltzmann
        if params.simulationType == "Langevin":
            
            # Kinetic energy figure
            plt.figure(kineticEnergyFig.number)
            
            # Get the Thermal equilibration limits
            start = -1
            end = -1
            for stage in meas.stage:
                if stage["Name"] == "ThermalLJEquilibration":
                    start = stage["Frame"][0]
                    end = stage["Frame"][1]
            
            if start != -1:
                
                # Plot measured
                plt.plot(meas.time[start:end], meas.kineticEnergy[start:end], "-bx", label="Kinetic Energy")

                # Plot expected (half kT per degree of freedom)
                xVals = np.linspace(0,plt.xlim()[1],100)
                yVals = 1.5 * params.kbt * params.numParticles * np.ones(100, dtype="float")
                plt.plot(xVals, yVals, "--k", label="Equipartition")

                # Prettify
                plt.xlabel("Time (ns)", fontsize=MEDIUM)
                plt.ylabel("Kinetic Energy (pN.nm)", fontsize=MEDIUM)
                plt.legend(loc=0, fontsize=SMALL)

                oFname = os.path.splitext(params.mFname_Equilibration)[0] + "_ThermalLJEquilKineticEnergy.png"
                plt.savefig(oFname, dpi=100)


        # LJ energy
        plt.figure(ljEnergyFig.number)
        
        # Get the equilibration limits
        start = -1
        end = -1
        for stage in meas.stage:
            if stage["Name"] == "ThermalLJEquilibration":
                start = stage["Frame"][0]
                end = stage["Frame"][1]

        if start != -1:

            # Plot measured
            plt.plot(meas.time[start:end], meas.ljEnergy[start:end], "-bx", label="LJ Energy")

            # Expected unknown, but should be stable
 
            # Prettify
            plt.xlabel("Time (ns)", fontsize=MEDIUM)
            plt.ylabel("LJ Energy (pN.nm)", fontsize=MEDIUM)
            plt.legend(loc=0, fontsize=SMALL)

            oFname = os.path.splitext(params.mFname_Equilibration)[0] + "_ThermalLJEquilLJEnergy.png"
            plt.savefig(oFname, dpi=100)           

        plt.show()
        plt.clf()

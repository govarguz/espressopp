#!/usr/bin/env python                                                               
# -*- coding: iso-8859-1 -*-                                                        

###########################################################################
#                                                                         #
#  ESPResSo++ Benchmark Python script for a Lennard Jones System          #
#                                                                         #
###########################################################################

import time
import espressopp
import csv

nsteps      = 10
isteps      = 10
rc          = pow(2.0, 1.0/6.0)
skin        = 0.4
timestep    = 0.005

# set temperature to None for NVE-simulations
temperature = 1.0

######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################
print espressopp.Version().info()
print 'Setting up simulation ...'
bonds, angles, x, y, z, Lx, Ly, Lz = espressopp.tools.convert.lammps.read('polymer_melt.lammps')
bonds, angles, x, y, z, Lx, Ly, Lz = espressopp.tools.replicate(bonds, angles, x, y, z, Lx, Ly, Lz, xdim=2, ydim=1, zdim=1)
num_particles = len(x)
density = num_particles / (Lx * Ly * Lz)
box = (Lx, Ly, Lz)
system, integrator = espressopp.standard_system.Default(box=box, rc=rc, skin=skin, dt=timestep, temperature=temperature)

# add particles to the system and then decompose
# do this in chunks of 1000 particles to speed it up
props = ['id', 'type', 'mass', 'pos']
new_particles = []
for i in range(num_particles):
  part = [i + 1, 0, 1.0, espressopp.Real3D(x[i], y[i], z[i])]
  new_particles.append(part)
  if i % 1000 == 0:
    system.storage.addParticles(new_particles, *props)
    system.storage.decompose()
    new_particles = []
system.storage.addParticles(new_particles, *props)

# HHack: Added NPartPerMPIrank
partNpMPIcp1=espressopp.analysis.NRealPart(system).compute()
#print espressopp.analysis.NRealPart(system).compute()
system.storage.decompose()

# Lennard-Jones with Verlet list
vl      = espressopp.VerletList(system, cutoff = rc)
potLJ   = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=0)
interLJ = espressopp.interaction.VerletListLennardJones(vl)
interLJ.setPotential(type1=0, type2=0, potential=potLJ)
system.addInteraction(interLJ)

# FENE bonds
fpl = espressopp.FixedPairList(system.storage)
fpl.addBonds(bonds)
potFENE = espressopp.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
interFENE = espressopp.interaction.FixedPairListFENE(system, fpl, potFENE)
system.addInteraction(interFENE)

# Cosine with FixedTriple list
ftl = espressopp.FixedTripleList(system.storage)
ftl.addTriples(angles)
potCosine = espressopp.interaction.Cosine(K=1.5, theta0=3.1415926)
interCosine = espressopp.interaction.FixedTripleListCosine(system, ftl, potCosine)
system.addInteraction(interCosine)

# print simulation parameters
print ''
print 'number of particles = ', num_particles
print 'density             = ', density
print 'rc                  = ', rc
print 'dt                  = ', integrator.dt
print 'skin                = ', system.skin
print 'temperature         = ', temperature
print 'nsteps              = ', nsteps
print 'isteps              = ', isteps
print 'NodeGrid            = ', system.storage.getNodeGrid()
print 'CellGrid            = ', system.storage.getCellGrid()
print ''

# espressopp.tools.decomp.tuneSkin(system, integrator)
'''
Saving Number of Interactions per core
Filename: "NIntPerRank.dlb"
example with 4 cores
csv llike:118747 118747 118477 118477
Filename: "NPartsPerRank.dlb"
csv llike:118747 118747 118477 118477

'''
		
f1=open("NIntPerRank.dlb","wt")		
wr1=csv.writer(f1,delimiter=" ")
pt1=[]
f2=open("NPartsPerRank.dlb","wt")		
wr2=csv.writer(f2,delimiter=" ")
pt2=[]
f3=open("PartsHangoverPerRankPerDir.dlb","wt")		
wr3=csv.writer(f3,delimiter=" ")
pt3=[]
f4=open("PartsCommCellsPerRankPerDir.dlb","wt")		
wr4=csv.writer(f4,delimiter=" ")
pt4=[]
'''
Finished creating load balancing ranks
'''
system.storage.resetParticlesHangoverCounters()
system.storage.resetParticlesCommCellsCounters()
espressopp.tools.analyse.info(system, integrator)
start_time = time.clock()
for k in range(nsteps):
  integrator.run(isteps)
  l=len(system.getInteraction(0).getVerletList().getAllPairs())
  p=[]
  for k2 in range(l):    	
    p.append(len(system.getInteraction(0).getVerletList().getAllPairs()[k2]))
  pt1.append(p)
  pt2.append(espressopp.analysis.NRealPart(system).compute())
  pt3.append(system.storage.getCounters())
  pt4.append(system.storage.getCommCellsTotalParticles())
  system.storage.resetParticlesHangoverCounters()
  system.storage.resetParticlesCommCellsCounters()
  espressopp.tools.analyse.info(system, integrator)
wr1.writerows(pt1)
wr2.writerows(pt2)
wr3.writerows(pt3)
wr4.writerows(pt4)
f1.close()
f2.close()
f3.close()
f4.close()
end_time = time.clock()
espressopp.tools.analyse.info(system, integrator)
espressopp.tools.analyse.final_info(system, integrator, vl, start_time, end_time)
'''
Saving Number of Interactions per core
Filename: "NIntPerRank.dlb"
csv llike: 
'''
#partNpMPIcp2=espressopp.analysis.NPart(system).compute()
#with open("NpartsPerRank.dlb","w") as f:
#  wr1=csv.writer(f,delimiter=" ")
#  wr1.writerows([partNpMPIcp1,partNpMPIcp2])
#fileout1.close()
#print espressopp.analysis.NRealPart(system).compute()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from crpropa import *


df=pd.read_csv('Cut5EeV.csv') # Defalut Energy Cut is 5 EeV
ID=df['AugerID'].values
E=df['E'].values
l=df['L'].values
b=df['B'].values

print(len(ID))

print(np.max(l),np.min(l))
print(np.max(b),np.min(b))
print(np.max(E))


#From here Example starts.

Field = JF12Field()
pid= - nucleusId(1,1) # pid should be inverse to make charge time symmetry.
# simulation setup
sim = ModuleList()
sim.add(PropagationCK(Field, 1e-8, 0.5 * parsec, 15 * parsec)) #Cash-Karp integrator.
obs = Observer()
obs.add(ObserverSurface( Sphere(Vector3d(0), 20 * kpc) ))
position = Vector3d(-8.5, 0, 0) * kpc
# obs.onDetection(TextOutput('galactic_backtracking.txt', Output.Event3D))
sim.add(obs)

# Additional parameters may be introduced depending on the direction of the simulation. (However, we did not introduce it in this analysis.)
"""
sim.add(Redshift())
sim.add(PhotoPionProduction(CMB()))
sim.add(PhotoPionProduction(IRB()))
sim.add(PhotoDisintegration(CMB()))
sim.add(PhotoDisintegration(IRB()))
sim.add(NclearDecay())
sim.add(ElectronPairProduction(CMB()))
sim.add(ElectronPairProduction(IRB()))
"""
wf=open('JF12_5EeVCut.csv','w')

head="AugerID,E,L,B,dL,dB\n"
wf.write(head)

#dl,db=[],[]
start=time.time()
for i in range(len(ID)):
    direction=Vector3d()
    lon=l[i]
    lat=np.pi/2-b[i]
    Energy=E[i]
    direction.setRThetaPhi(1,lat,lon)
    p=ParticleState(pid,Energy,position,direction)
    c=Candidate(p)
    sim.run(c)
    d1 = c.current.getDirection()
    #dl.append(d1.getPhi())
    #db.append(d1.getTheta())
    dl=d1.getPhi()
    db=np.pi/2-d1.getTheta()
    line="%d,%lf,%lf,%lf,%lf,%lf\n"%(ID[i],E[i],l[i],b[i],dl,db)
    wf.write(line)
end=time.time()

print(end-start)
wf.close()

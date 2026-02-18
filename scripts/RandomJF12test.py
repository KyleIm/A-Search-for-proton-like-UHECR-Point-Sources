import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from crpropa import *

E=5
l=0
b=0

#l=0
#b=0

dl_l=[]
db_l=[]


for i in range(1000):
    Field = JF12Field(True)
    pid= - nucleusId(1,1) # pid should be inverse to make charge time symmetry.
    sim = ModuleList()
    sim.add(PropagationCK(Field, 1e-8, 0.5 * parsec, 15 * parsec)) #Cash-Karp integrator.
    obs = Observer()
    obs.add(ObserverSurface( Sphere(Vector3d(0), 20 * kpc) ))
    direction=Vector3d()
    lon=l
    lat=np.pi/2-b
    Energy=E
    position = Vector3d(-8.5, 0, 0) * kpc
    sim.add(obs)
    direction.setRThetaPhi(1,lat,lon)
    p=ParticleState(pid,Energy,position,direction)
    c=Candidate(p)
    sim.run(c)
    d1 = c.current.getDirection()
    dl=d1.getPhi()
    db=np.pi/2-d1.getTheta()
    dl_l.append(dl)
    db_l.append(db)

Field = JF12Field(True)
pid= - nucleusId(1,1) # pid should be inverse to make charge time symmetry.
sim = ModuleList()
sim.add(PropagationCK(Field, 1e-8, 0.5 * parsec, 15 * parsec)) #Cash-Karp integrator.
obs = Observer()
obs.add(ObserverSurface( Sphere(Vector3d(0), 20 * kpc) ))
direction=Vector3d()
lon=l
lat=np.pi/2-b
Energy=E
position = Vector3d(-8.5, 0, 0) * kpc
sim.add(obs)
direction.setRThetaPhi(1,lat,lon)
p=ParticleState(pid,Energy,position,direction)
c=Candidate(p)
sim.run(c)
d1 = c.current.getDirection()
dl=d1.getPhi()
db=np.pi/2-d1.getTheta()


dL=np.array(dl_l)
dB=np.array(db_l)


plt.figure()

ax = plt.subplot(111, projection='mollweide')
ax.grid(True)

ax.scatter(dL*-1, dB, s=15, label='randfield')
ax.scatter(l*-1, b, s=40, c='g', label='Origin')
ax.scatter(dl*-1, db, s=40, c='r', label='JF12')

# Change axis tick font size
xticks = np.radians([-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180])
ax.set_xticks(xticks)
ax.set_xticklabels(["180°", "150°", "120°", "90°", "60°", "30°", "0°", "-30°", "-60°", "-90°", "-120°", "-150°", "-180°"], fontsize=20)

#ax.tick_params(axis='x', labelsize=20)  # X-axis tick font size
ax.tick_params(axis='y', labelsize=20)  # Y-axis tick font size
ax.set_xlabel("L", fontsize=25)
ax.set_ylabel("B", fontsize=25)
ax.legend(fontsize=25)
plt.title("Deflection (0, 0), 5EeV", fontsize=30)
plt.show()

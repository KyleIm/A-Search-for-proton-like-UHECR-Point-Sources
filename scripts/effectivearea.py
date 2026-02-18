import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from healpy.newvisufunc import projview
from crpropa import *
import time
import pandas as pd

def Angdistance(l1,l2,b1,b2): #Analytic formula for distance on sphere
    D=np.rad2deg(np.arccos(np.sin(b1)*np.sin(b2) + np.cos(b1)*np.cos(b2)*np.cos(l2-l1)))   
    return D

def dLdB(E,l,b,Jdl,Jdb,trial,D_l):
    for i in range(trial):
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
        # obs.onDetection(TextOutput('galactic_backtracking.txt', Output.Event3D))
        sim.add(obs)
        direction.setRThetaPhi(1,lat,lon)
        p=ParticleState(pid,Energy,position,direction)
        c=Candidate(p)
        sim.run(c)
        d1 = c.current.getDirection()
        dl=d1.getPhi()
        db=np.pi/2-d1.getTheta()
        D = Angdistance(dl,Jdl,db,Jdb)
        D_l.append(D)

start=time.time()


df=pd.read_csv('Cut50EeV.csv')
l=df['L'].values
b=df['B'].values

E = 5

Distance_l = []

for i in range(len(l)):
    Field = JF12Field() # Default false
    pid= - nucleusId(1,1) # pid should be inverse to make charge time symmetry.
    sim = ModuleList()
    sim.add(PropagationCK(Field, 1e-8, 0.5 * parsec, 15 * parsec)) #Cash-Karp integrator.
    obs = Observer()
    obs.add(ObserverSurface( Sphere(Vector3d(0), 20 * kpc) ))
    direction=Vector3d()
    lon=l[i]
    lat=np.pi/2-b[i]
    Energy=E
    position = Vector3d(-8.5, 0, 0) * kpc
    # obs.onDetection(TextOutput('galactic_backtracking.txt', Output.Event3D))
    sim.add(obs)
    direction.setRThetaPhi(1,lat,lon)
    p=ParticleState(pid,Energy,position,direction)
    c=Candidate(p)
    sim.run(c)
    d1 = c.current.getDirection()
    Jdl=d1.getPhi()
    Jdb=np.pi/2-d1.getTheta()
    dLdB(E,l[i],b[i],Jdl,Jdb,1000,Distance_l)
            
end=time.time()
print(end-start)
print(np.mean(Distance_l),np.std(Distance_l))


mean_distance = np.mean(Distance_l)
std_distance = np.std(Distance_l)

# Plot histogram
plt.figure(figsize=(10, 6))
plt.hist(Distance_l, bins=30, edgecolor='blue', alpha=0.7, label="Distances")

# Add vertical lines for mean and ±1 sigma
plt.axvline(mean_distance, color='red', linestyle='--', linewidth=3, label=f'Mean : {mean_distance:.2f}')
#plt.axvline(mean_distance - std_distance, color='green', linestyle='--', linewidth=1.5, label='-1σ')
plt.axvline(mean_distance + std_distance, color='green', linestyle='--', linewidth=1.5, label=f'+1σ : {mean_distance + std_distance:.2f}')


# Customize plot
plt.title("Effective Distance for 5EeV", fontsize=30)
plt.xlabel("Angular Distance (°)", fontsize=25)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(fontsize=20)
plt.yscale('log')
#plt.grid(alpha=0.3)

# Show plot
plt.show()

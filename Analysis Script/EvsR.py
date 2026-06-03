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

ENERGY_MIN = 1.0
ENERGY_MAX = 10.0
ENERGY_STEP = 0.5
TRIALS_PER_EVENT = 10

def deterministic_seed(stream, energy_idx, event_idx, trial_idx):
    seed = (
        2166136261
        ^ (stream + 1) * 16777619
        ^ (energy_idx + 11) * 2246822519
        ^ (event_idx + 101) * 3266489917
        ^ (trial_idx + 1009) * 668265263
    ) & 0xFFFFFFFF
    return seed if seed != 0 else 1

def dLdB(E,l,b,Jdl,Jdb,trial,D_l,energy_idx,event_idx):
    max_attempts = 5
    for i in range(trial):
        attempt_count = 0
        while attempt_count < max_attempts:
            Random.instance().seed(deterministic_seed(1, energy_idx, event_idx, i * max_attempts + attempt_count))
            Field = JF12Field(True)
            pid= - nucleusId(1,1) # pid should be inverse to make charge time symmetry.
            sim = ModuleList()
            sim.add(PropagationCK(Field, 1e-8, 0.5 * parsec, 15 * parsec)) #Cash-Karp integrator.
            obs = Observer()
            obs.add(ObserverSurface( Sphere(Vector3d(0), 20 * kpc) ))
            sim.add(obs)
            
            direction=Vector3d()
            lon=l
            lat=np.pi/2-b
            Energy=E
            position = Vector3d(-8.5, 0, 0) * kpc
            
            direction.setRThetaPhi(1,lat,lon)
            p=ParticleState(pid,Energy,position,direction)
            c=Candidate(p)
            sim.run(c)
            
            d1 = c.current.getDirection()
            dl=d1.getPhi()
            db=np.pi/2-d1.getTheta()
            
            if np.isfinite(dl) and np.isfinite(db):
                D = Angdistance(dl,Jdl,db,Jdb)
                D_l.append(D)
                break
            else:
                attempt_count += 1
                print(attempt_count, np.nan)
        
        if attempt_count == max_attempts:
            D_l.append(np.nan)

start=time.time()

print(start)
df=pd.read_csv('Cut50EeV.csv')
l=df['L'].values
b=df['B'].values

print(len(l))

N_3 = np.rad2deg(np.arccos(1-1/(6*3**2)))
N_6 = np.rad2deg(np.arccos(1-1/(6*6**2)))
N_9 = np.rad2deg(np.arccos(1-1/(6*9**2)))
N_12 = np.rad2deg(np.arccos(1-1/(6*12**2)))
N_15 = np.rad2deg(np.arccos(1-1/(6*15**2)))
N_18 = np.rad2deg(np.arccos(1-1/(6*18**2)))

Mean_cut = []
E_l = []

Distance_l = []
energies = np.arange(ENERGY_MIN, ENERGY_MAX + 0.5 * ENERGY_STEP, ENERGY_STEP)
with open("EvsR.csv",'w') as f:
    head = "E, mean, sigma \n"
    f.write(head)
    for i, E in enumerate(energies):
        for j in range(len(l)):
            Random.instance().seed(deterministic_seed(0, i, j, 0))
            Field = JF12Field() # Default false
            pid= - nucleusId(1,1) # pid should be inverse to make charge time symmetry.
            sim = ModuleList()
            sim.add(PropagationCK(Field, 1e-8, 0.5 * parsec, 15 * parsec)) #Cash-Karp integrator.
            obs = Observer()
            obs.add(ObserverSurface( Sphere(Vector3d(0), 20 * kpc) ))
            direction=Vector3d()
            lon=l[j]
            lat=np.pi/2-b[j]
            #Energy=E
            position = Vector3d(-8.5, 0, 0) * kpc
            # obs.onDetection(TextOutput('galactic_backtracking.txt', Output.Event3D))
            sim.add(obs)
            direction.setRThetaPhi(1,lat,lon)
            p=ParticleState(pid,E,position,direction)
            c=Candidate(p)
            sim.run(c)
            d1 = c.current.getDirection()
            Jdl=d1.getPhi()
            Jdb=np.pi/2-d1.getTheta()
            dLdB(E,l[j],b[j],Jdl,Jdb,TRIALS_PER_EVENT,Distance_l,i,j)
        Mean_cut.append(np.mean(Distance_l))
        print(np.mean(Distance_l),np.mean(Distance_l))
        E_l.append(E)
        f.write("%lf, %lf, %lf \n" %(E, np.mean(Distance_l), np.std(Distance_l)))
        Distance_l = []
    
            
end=time.time()
print(end-start)


# Plot histogram
plt.figure(figsize=(10, 6))
plt.scatter(E_l,Mean_cut,c='r',s=30,label='mean distance')

plt.axhline(N_3, color='orange', linestyle='--', label=r"$N_{\mathrm{side}}$=3")
plt.axhline(N_6, color='yellow', linestyle='--', label=r"$N_{\mathrm{side}}$=6")
plt.axhline(N_9, color='green', linestyle='--', label=r"$N_{\mathrm{side}}$=9")
plt.axhline(N_12, color='indigo', linestyle='--', label=r"$N_{\mathrm{side}}$=12")
plt.axhline(N_15, color='magenta', linestyle='--', label=r"$N_{\mathrm{side}}$=15")
plt.axhline(N_18, color='purple', linestyle='--', label=r"$N_{\mathrm{side}}$=18")


# Customize plot
plt.title("Energy vs Radius Cut", fontsize=40)
plt.xlabel("Energy (EeV)",fontsize=35)
plt.ylabel("Radius (°)", fontsize=35)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.legend(fontsize=25)
#plt.grid(alpha=0.3)

# Show plot
plt.show()

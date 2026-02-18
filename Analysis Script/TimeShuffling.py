import numpy as np
import pandas as pd
from datetime import datetime
from crpropa import *
import time as ti

J2000=datetime.utcfromtimestamp(946684800)
lat_0=np.deg2rad(-35.2) #South
lon_0=np.deg2rad(-69.3) #West
Dec_G=np.deg2rad(27.12825)
RA_G=np.deg2rad(192.85948)
l_NCP=np.deg2rad(122.93192)

def Hor2Equ(azi,zen,t): # phi, theta, t (time in UTC)
    azin=np.pi/2-azi
    Dec=np.arcsin(np.cos(zen)*np.sin(lat_0)+np.sin(zen)*np.cos(lat_0)*np.cos(azin))
    nume=-1*np.sin(azin)*np.sin(zen) #y
    deno=np.cos(zen)*np.cos(lat_0)-np.sin(zen)*np.sin(lat_0)*np.cos(azin) #x
    H=np.arctan2(nume,deno) #tan(y/x)=np.arctan2(y,x)
    if(t-H<0):
        RA=t-H+2*np.pi
    elif(t-H>2*np.pi):
        RA=t-H-2*np.pi
    else:
        RA=t-H
    return(RA,Dec)

def Equ2Gal(RA,Dec):
    b=np.arcsin(np.sin(Dec)*np.sin(Dec_G)+np.cos(Dec)*np.cos(Dec_G)*np.cos(RA-RA_G))
    deno=np.sin(Dec)*np.cos(Dec_G)-np.cos(Dec)*np.sin(Dec_G)*np.cos(RA-RA_G) #x
    nume=np.cos(Dec)*np.sin(RA-RA_G) #y
    l=l_NCP-np.arctan2(nume,deno)
    if(l>np.pi):
        l=l-2*np.pi
    return (l,b)

def UTC2LST(UTC): #Time should be given as epoch
    YMDHMS=datetime.utcfromtimestamp(UTC)
    strUTC=str(YMDHMS)
    UTH=float(strUTC[-8:-6])
    UTM=float(strUTC[-5:-3])
    UTS=float(strUTC[-2:])
    UT=UTH+UTM/60+UTS/3600
    d=(YMDHMS-J2000).days
    val=100.46+0.985647*d-69.3+15*UT 
    LSTdeg=val%360 #val is in deg but we need leftover of large deg.
    H=LSTdeg//15
    M=(LSTdeg-H*15)//0.25
    S=(LSTdeg-H*15-M*0.25)//(1/240)
    LSTrad=np.deg2rad(LSTdeg)
    return LSTrad


def Deflect(l,b,E):
    check = 0
    while(check == 0):
        Field = JF12Field()
        pid= - nucleusId(1,1) # pid should be inverse to make charge time symmetry.
        # simulation setup
        sim = ModuleList()
        sim.add(PropagationCK(Field, 1e-4, 0.1 * parsec, 100 * parsec)) #Cash-Karp integrator.
        obs = Observer()
        obs.add(ObserverSurface(Sphere(Vector3d(0), 20 * kpc) ))
        position = Vector3d(-8.5, 0, 0) * kpc
        sim.add(obs)
        direction=Vector3d()
        lon=l
        lat=np.pi/2-b
        Energy=E
        direction.setRThetaPhi(1,lat,lon)
        p=ParticleState(pid,Energy,position,direction)
        c=Candidate(p)
        sim.run(c)
        d1 = c.current.getDirection()
        dl=d1.getPhi()
        db=np.pi/2-d1.getTheta()
        if np.isfinite(dl) and np.isfinite(db):
            check = 1

    return (dl,db)

def Shuffle(t):
    ts=t.copy()
    np.random.shuffle(ts)
    return (ts)

df=pd.read_csv('Cut5EeV.csv')
ID=df['AugerID'].values
Theta=df['Theta'].values
Phi=df['Phi'].values
RA=df['RA'].values
Dec=df['Dec'].values
l=df['L'].values
b=df['B'].values
UTC=df['UTC'].values
E=df['E'].values

print(np.max(E))


f=open('test5EeV.csv','w')
head="Trial,RA,Dec,L,B,dL,dB,E\n"
f.write(head)

start=ti.time()
for i in range(20):
    sUTC=Shuffle(UTC)
    for j in range(len(ID)):
        time=UTC2LST(sUTC[j])
        ra,dec=Hor2Equ(Phi[j],Theta[j],time)
        L,B=Equ2Gal(ra,dec)
        dl,db=Deflect(L,B,E[j])
        Trial=i*len(ID)+j
        data="%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n"%(Trial,ra,dec,L,B,dl,db,E[j])
        f.write(data)
        
end=ti.time()

print(end-start)
print(len(ID))
f.close()

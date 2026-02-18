import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
from healpy.newvisufunc import projview
import time as ti

df1=pd.read_csv('Cut5EeV.csv')
ID=df1['AugerID'].values
L=df1['L'].values
B=np.pi/2 - df1['B'].values

df=pd.read_csv('Background_psf.csv')
ucback = df['ucback_psf_map'].values
back=df['back_psf_map'].values


df3=pd.read_csv('JF12_5EeVCut.csv')
ID3=df3['AugerID'].values
dL=df3['dL'].values
dB=np.pi/2 - df3['dB'].values

psf_map = back/np.sum(back) *(43463*20)

def LiMa(n,b,alpha):
    if(b==n):
        lima=np.sqrt(2*(n*np.log((n+alpha*n)/(b+alpha*n))+b/alpha*np.log((b+alpha*b)/(b+alpha*n))))
        print("Something went wrong")
    elif(n==0):
        lima=(n-b)/abs(n-b)*np.sqrt(2*(b/alpha*np.log((b+alpha*b)/(b+alpha*n)))) #Not making inside ln tobe 0)
    else:
        lima=(n-b)/abs(n-b)*np.sqrt(2*(n*np.log((n+alpha*n)/(b+alpha*n))+b/alpha*np.log((b+alpha*b)/(b+alpha*n))))
    if(n*np.log((n+alpha*n)/(b+alpha*n))+b/alpha*np.log((b+alpha*b)/(b+alpha*n))<0):
        print("Imaginary")
    return lima


def Angdistance(l1,l2,b1,b2): #Analytic formula for distance on sphere
    D=np.arccos(np.cos(b1)*np.cos(b2) + np.sin(b1)*np.sin(b2)*np.cos(l2-l1))
    return D

def PSF(D, sigma):
    psf = np.exp(-1*D**2/(2*sigma**2))/(2*np.pi*sigma**2)
    return psf
    
nside=128
sigma = np.deg2rad(2.16)
npix=hp.nside2npix(nside)

start = ti.time()

theta, phi = hp.pix2ang(nside, np.arange(npix))

ucorg_psf_map = np.zeros(npix, dtype=float)
org_psf_map = np.zeros(npix, dtype=float)


print(ti.time() - start)

dglimamap_l=[]

f = open('Realpsf.csv','w')
head = "pixel,uc_psf_map,back_psf_map\n"
f.write(head)
pixel = 0

for i in range(len(dL)):
    D = Angdistance(L[i], phi, B[i], theta)
    ucorg_psf_map += PSF(D, sigma)

print("UC complete")
    
for j in range(len(dL)):
    D = Angdistance(dL[j], phi, dB[j], theta)
    org_psf_map += PSF(D, sigma)

print("C complete")
    
for k in range(len(org_psf_map)):
    line = "%d,%lf,%lf\n"%(pixel,ucorg_psf_map[k],org_psf_map[k])
    f.write(line)
    pixel += 1

f.close()

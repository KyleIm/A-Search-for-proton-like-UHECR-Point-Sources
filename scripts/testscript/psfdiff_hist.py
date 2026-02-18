import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
from healpy.newvisufunc import projview


df=pd.read_csv('../JF12_5EeVCut.csv')
l=df['L'].values
b=np.pi/2-df['B'].values
dl=df['dL'].values
db=np.pi/2-df['dB'].values

df2=pd.read_csv('../test5EeV.csv')
bl=df2['L'].values
bb=np.pi/2-df2['B'].values
bdl=df2['dL'].values
bdb=np.pi/2-df2['dB'].values



def Angdistance(l1,l2,b1,b2): #Analytic formula for distance on sphere
    D=np.arccos(np.cos(b1)*np.cos(b2) + np.sin(b1)*np.sin(b2)*np.cos(l2-l1))
    return D

def PSF(D, sigma):
    psf = np.exp(-1*D**2/(2*sigma**2))/(2*np.pi*sigma**2)
    return psf

nside=15
sigma = np.deg2rad(2.16)
npix=hp.nside2npix(nside) #12*nside**2

theta, phi = hp.pix2ang(nside, np.arange(npix))

org_psf_map = np.zeros(npix, dtype=float)
jf12_psf_map = np.zeros(npix, dtype=float)

back_psf_map = np.zeros(npix, dtype=float)
backjf12_psf_map = np.zeros(npix, dtype=float)

for i in range(len(l)):
    D = Angdistance(l[i], phi, b[i], theta)
    org_psf_map += PSF(D, sigma)
    dD = Angdistance(dl[i], phi, db[i], theta)
    jf12_psf_map += PSF(dD, sigma)

for j in range(len(bl)):
    bD = Angdistance(bl[j], phi, bb[j], theta)
    back_psf_map += PSF(bD, sigma)
    bdD = Angdistance(bdl[j], phi, bdb[j], theta)
    backjf12_psf_map += PSF(bdD, sigma)


org_mask = org_psf_map/np.sum(org_psf_map) * 43463 > 0.5
hpmap = (org_psf_map/np.sum(org_psf_map) - back_psf_map/np.sum(back_psf_map)) * 43463
jf12_mask = jf12_psf_map/np.sum(jf12_psf_map) * 43463 > 0.5
bhpmap = (jf12_psf_map/np.sum(jf12_psf_map) - backjf12_psf_map/np.sum(backjf12_psf_map)) * 43463



f_hpmap = hpmap[org_mask]
bf_hpmap = bhpmap[jf12_mask]

MIN = min(np.min(f_hpmap),np.min(bf_hpmap))
MAX = max(np.max(f_hpmap),np.max(bf_hpmap))

plt.hist(bf_hpmap,bins=30,range=(MIN - 0.5,MAX + 0.5),facecolor='None',lw=3,edgecolor='blue',label='JF12 : %d'%len(bf_hpmap))
plt.hist(f_hpmap,bins=30,range=(MIN - 0.5,MAX + 0.5),facecolor='None',lw=2,edgecolor='red',label='RAW : %d'%len(f_hpmap))
plt.legend()
plt.title("Histogram for PSF Difference, E > 5EeV, nside = 15, - log",fontsize=30)
plt.yscale('log')
plt.show()

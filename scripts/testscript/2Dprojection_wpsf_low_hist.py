import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
from healpy.newvisufunc import projview


df=pd.read_csv('../test5EeV.csv')
l=df['L'].values
b=np.pi/2-df['B'].values
dl=df['dL'].values
db=np.pi/2-df['dB'].values


def Angdistance(l1,l2,b1,b2): #Analytic formula for distance on sphere
    D=np.arccos(np.cos(b1)*np.cos(b2) + np.sin(b1)*np.sin(b2)*np.cos(l2-l1))
    return D

def PSF(D, sigma):
    psf = D/(sigma**2)*np.exp(-1*D**2/(2*sigma**2))/(2*np.pi*sigma**2) / (2*np.pi*np.sin(D))
    return psf

nside=15
sigma = np.deg2rad(2.16)
npix=hp.nside2npix(nside) #12*nside**2

theta, phi = hp.pix2ang(nside, np.arange(npix))

org_psf_map = np.zeros(npix, dtype=float)
jf12_psf_map = np.zeros(npix, dtype=float)


for i in range(len(l)):
    D = Angdistance(l[i], phi, b[i], theta)
    org_psf_map += PSF(D, sigma)
    dD = Angdistance(dl[i], phi, db[i], theta)
    jf12_psf_map += PSF(dD, sigma)

hpmap = org_psf_map/np.sum(org_psf_map) * 43463 * 20
dhpmap = jf12_psf_map/np.sum(jf12_psf_map) * 43463 * 20

f_hpmap = hpmap[hpmap > 0.5]
df_hpmap = dhpmap[dhpmap > 0.5]

MIN = min(np.min(f_hpmap),np.min(df_hpmap))
MAX = max(np.max(f_hpmap),np.max(df_hpmap))

plt.hist(df_hpmap,bins=30,range=(MIN - 0.5,MAX + 0.5),facecolor='None',lw=3,edgecolor='blue',label='JF12 : %d'%len(df_hpmap))
plt.hist(f_hpmap,bins=30,range=(MIN - 0.5,MAX + 0.5),facecolor='None',lw=2,edgecolor='red',label='RAW : %d'%len(f_hpmap))
plt.legend()
plt.title("Histogram for PSF Backgrounds, E > 5EeV, nside = 15, - log",fontsize=30)
plt.yscale('log')
plt.show()

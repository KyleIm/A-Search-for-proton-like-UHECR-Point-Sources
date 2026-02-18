import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
from healpy.newvisufunc import projview


df=pd.read_csv('JF12_5EeVCut.csv')
l=df['L'].values
b=np.pi/2-df['B'].values

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
    #simple=(sig-background)/np.sqrt(background)
    #return (simple)

def Angdistance(l1,l2,b1,b2): #Analytic formula for distance on sphere
    D=np.arccos(np.cos(b1)*np.cos(b2) + np.sin(b1)*np.sin(b2)*np.cos(l2-l1))
    return D

def PSF(D, sigma):
    psf = D/(sigma**2)*np.exp(-1*D**2/(2*sigma**2))/(2*np.pi*sigma**2)/(2*np.pi*np.sin(D))
    return psf

nside=15
sigma = np.deg2rad(2.16)
npix=hp.nside2npix(nside) #12*nside**2

theta, phi = hp.pix2ang(nside, np.arange(npix))

org_psf_map = np.zeros(npix, dtype=float)

for i in range(len(l)):
    D = Angdistance(l[i], phi, b[i], theta)
    org_psf_map += PSF(D, sigma)

print(PSF(D, sigma))
    
hpmap = org_psf_map/np.sum(org_psf_map) * 43463

projview(
    hpmap,
    coord=["G"],
    graticule=True,
    graticule_labels=True,
    #title="Binning with healpix",
    xlabel="l",
    ylabel="b",
    flip="astro",
    projection_type="mollweide",
)

ax = plt.gca()
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(20)  

ax.set_xlabel("L", fontsize=25)
ax.set_ylabel("B", fontsize=25)
    
cbar = plt.gcf().axes[-1]
cbar.tick_params(labelsize=20)
    

plt.title("2D projection of PSF (RAW), E < 5EeV, nside = 15",fontsize=30)
plt.show()

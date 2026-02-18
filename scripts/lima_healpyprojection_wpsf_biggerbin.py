import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
from healpy.newvisufunc import projview

df1=pd.read_csv('Cut5EeV.csv')
ID=df1['AugerID'].values
L=df1['L'].values
B=np.pi/2 - df1['B'].values

df2=pd.read_csv('test5EeV.csv')
Trial=df2['Trial'].values
bL=df2['L'].values
bB=np.pi/2 - df2['B'].values
bdL=df2['dL'].values
bdB=np.pi/2 - df2['dB'].values

df3=pd.read_csv('JF12_5EeVCut.csv')
ID3=df3['AugerID'].values
dL=df3['dL'].values
dB=np.pi/2 - df3['dB'].values

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
    psf = D * np.exp(-1*D**2/(2*sigma**2))/(2*np.pi*sigma**2)
    return psf

nside=15
sigma = np.deg2rad(2.16)
npix=hp.nside2npix(nside) #12*nside**2

theta, phi = hp.pix2ang(nside, np.arange(npix))

org_psf_map = np.zeros(npix, dtype=float)
back_psf_map = np.zeros(npix, dtype=float)

ucorg_psf_map = np.zeros(npix, dtype=float)
ucback_psf_map = np.zeros(npix, dtype=float)

for i in range(len(dL)):
    D = Angdistance(dL[i], phi, dB[i], theta)
    org_psf_map += PSF(D, sigma)

org_psf_map = org_psf_map/np.sum(org_psf_map) * 43463
print(org_psf_map)

for j in range(len(bdL)):
    D = Angdistance(bdL[j], phi, bdB[j], theta)
    back_psf_map += PSF(D, sigma)

back_psf_map = back_psf_map/np.sum(back_psf_map) * 43463 * 20    
lima_map = []


print(back_psf_map)

for k in range(len(org_psf_map)):
    if(back_psf_map[k]<1e-2):
        lima_map.append(np.nan)
    else:
        lima=LiMa(org_psf_map[k],back_psf_map[k]/20,1/20)
        lima_map.append(lima)

hpmap = np.array(lima_map)
        
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

#print(hpmap)
a=np.where(hpmap==np.max(hpmap))
print(a)
print(hpmap[a])
print(len(hpmap))
#print(hpmap,np.max(hpmap))
plt.title("LiMa projection (JF12) with PSF : bins=2700, alpha=1/20,  E > 5EeV",fontsize=30)
#plt.figure()
#plt.subplot(111,projection='aitoff')
#plt.hexbin(dl,db)
#plt.hist(hpmap,bins=14)
#plt.yscale('log')
plt.show()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
from healpy.newvisufunc import projview

df=pd.read_csv('JF12_5EeVCut.csv')
l=df['dL'].values
b=np.pi/2-df['dB'].values

print(np.max(l),np.min(l),np.max(b),np.min(b))

nside=15
npix=hp.nside2npix(nside)

indices=hp.ang2pix(nside,b,l)

hpmap=np.zeros(npix)

for i in range(len(l)):
    hpmap[indices[i]]=hpmap[indices[i]]+1

count = 0
for j in range(len(hpmap)):
    (br,lr) = hp.pix2ang(nside,j)
    bd = 90 - np.rad2deg(br)
    ld = np.rad2deg(lr)
    if (bd < -15 and bd > -30 and hpmap[j] > 40):
        print(ld, bd, hpmap[j])
        count += 1

print(count)


projview(
    hpmap,
    coord=["G"],
    graticule=True,
    graticule_labels=True,
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
#plt.title("2D projection of Pierre Auger Raw Dataset, E > 5EeV, nside = 15",fontsize=30)
plt.title("2D projection of Pierre Auger JF12 Backtracked Dataset, E > 5EeV, nside = 15",fontsize=30)
#plt.title("2D projection of Time Shuffled Background (JF12 Backtracked), E > 5EeV, nside = 15",fontsize=30)


plt.show()

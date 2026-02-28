import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
from healpy.newvisufunc import projview

df=pd.read_csv('JF12_5EeVCut.csv')
l=df['L'].values
b=np.pi/2-df['B'].values


nside=15
npix=hp.nside2npix(nside)

indices=hp.ang2pix(nside,b,l)

hpmap=np.zeros(npix)

for i in range(len(l)):
    hpmap[indices[i]]=hpmap[indices[i]]+1

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
plt.title(r"2D projection of Pierre Auger Raw Dataset, E > 5EeV, $N_{side}$ = 15",fontsize=30)
#plt.title(r"2D projection of Pierre Auger JF12 Backtracked Dataset, E > 5EeV, $N_{side}$ = 15",fontsize=30)
#plt.title(r"2D projection of Time Shuffled Background Simulation, E > 5EeV, $N_{side}$ = 15", fontsize=30)



plt.show()

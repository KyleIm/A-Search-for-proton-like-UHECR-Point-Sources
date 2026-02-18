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


nside=15
npix=hp.nside2npix(nside)

indices=hp.ang2pix(nside,b,l)
dindices = hp.ang2pix(nside,db,dl)

hpmap = np.zeros(npix)
dhpmap = np.zeros(npix)

for i in range(len(l)):
    hpmap[indices[i]]=hpmap[indices[i]]+1

for j in range(len(dl)):
    dhpmap[dindices[j]]=dhpmap[dindices[j]]+1

f_hpmap = hpmap[hpmap > 0]
df_hpmap = dhpmap[dhpmap > 0]

MIN = min(np.min(f_hpmap),np.min(df_hpmap))
MAX = max(np.max(f_hpmap),np.max(df_hpmap))

plt.hist(df_hpmap,bins=30,range=(MIN - 0.5,MAX + 0.5),facecolor='None',lw=3,edgecolor='blue',label='JF12 : %d'%len(df_hpmap))
plt.hist(f_hpmap,bins=30,range=(MIN - 0.5,MAX + 0.5),facecolor='None',lw=2,edgecolor='red',label='RAW : %d'%len(f_hpmap))
plt.legend()
plt.title("Histogram for Background, E > 5EeV, nside = 15, - log",fontsize=30)
plt.yscale('log')
plt.show()

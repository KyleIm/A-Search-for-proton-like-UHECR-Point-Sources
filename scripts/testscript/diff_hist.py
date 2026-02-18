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

nside=15
npix=hp.nside2npix(nside)

indices=hp.ang2pix(nside,b,l)
dindices = hp.ang2pix(nside,db,dl)

bindices=hp.ang2pix(nside,bb,bl)
bdindices = hp.ang2pix(nside,bdb,bdl)

orgmap = np.zeros(npix)
jf12map = np.zeros(npix)

backmap = np.zeros(npix)
backjf12map = np.zeros(npix)



for i in range(len(l)):
    orgmap[indices[i]]=orgmap[indices[i]]+1
    jf12map[dindices[i]]=jf12map[dindices[i]]+1

for j in range(len(bl)):
    backmap[bindices[j]]=backmap[bindices[j]]+1
    backjf12map[bdindices[j]]=backjf12map[bdindices[j]]+1
    
hpfilter = orgmap > 0
hpmap = orgmap - backmap/20

dhpfilter = jf12map > 0
dhpmap = jf12map - backjf12map/20

print(dhpmap)

f_hpmap = hpmap[hpfilter]
df_hpmap = dhpmap[dhpfilter]

MIN = min(np.min(f_hpmap),np.min(df_hpmap))
MAX = max(np.max(f_hpmap),np.max(df_hpmap))

plt.hist(df_hpmap,bins=30,range=(MIN - 0.5,MAX + 0.5),facecolor='None',lw=3,edgecolor='blue',label='JF12 : %d'%len(df_hpmap))
plt.hist(f_hpmap,bins=30,range=(MIN - 0.5,MAX + 0.5),facecolor='None',lw=2,edgecolor='red',label='RAW : %d'%len(f_hpmap))
plt.legend()
plt.title("Histogram for Difference, E > 5EeV, nside = 15, - log",fontsize=30)
plt.yscale('log')
plt.show()

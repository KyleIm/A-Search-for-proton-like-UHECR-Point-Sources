import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
from healpy.newvisufunc import projview

df1=pd.read_csv('Cut5EeV.csv')
L_cut=df1['L'].values
B_cut=df1['B'].values

df1_syn=pd.read_csv('Sy_events_5EeV.csv')
L_syn=df1_syn['L'].values
B_syn=df1_syn['B'].values
sdL=df1_syn['dL'].values
sdB=df1_syn['dB'].values

L=np.concatenate((L_cut, L_syn))
B=np.concatenate((B_cut, B_syn))

df2=pd.read_csv('test5EeV.csv')
bL_test=df2['L'].values
bB_test=df2['B'].values
bdL_test=df2['dL'].values
bdB_test=df2['dB'].values

df2_syn=pd.read_csv('testSynthetic_5EeV.csv')
bL_syn=df2_syn['L'].values
bB_syn=df2_syn['B'].values
bdL_syn=df2_syn['dL'].values
bdB_syn=df2_syn['dB'].values

bL=np.concatenate((bL_test, bL_syn))
bB=np.concatenate((bB_test, bB_syn))
bdL=np.concatenate((bdL_test, bdL_syn))
bdB=np.concatenate((bdB_test, bdB_syn))

try:
    df3=pd.read_csv('JF12_5EeV.csv')
except FileNotFoundError:
    df3=pd.read_csv('JF12_5EeVCut.csv')
dL_cut=df3['dL'].values
dB_cut=df3['dB'].values

dL=np.concatenate((dL_cut, sdL))
dB=np.concatenate((dB_cut, sdB))

def LiMa(n,b,alpha):
    if(b==n):
        lima=np.sqrt(2*(n*np.log((n+alpha*n)/(b+alpha*n))+b/alpha*np.log((b+alpha*b)/(b+alpha*n))))
    elif(n==0):
        lima=(n-b)/abs(n-b)*np.sqrt(2*(b/alpha*np.log((b+alpha*b)/(b+alpha*n)))) #Not making inside ln tobe 0)
    else:
        lima=(n-b)/abs(n-b)*np.sqrt(2*(n*np.log((n+alpha*n)/(b+alpha*n))+b/alpha*np.log((b+alpha*b)/(b+alpha*n))))
    if(n*np.log((n+alpha*n)/(b+alpha*n))+b/alpha*np.log((b+alpha*b)/(b+alpha*n))<0):
        print("Imagninary number is returned") # When you see this message, it means something went wrong and at least one of the bin will not have Real Number result!
    return lima

nside=15
npix=hp.nside2npix(nside) #12*nside**2
ghpmap=np.zeros(npix)
dghpmap=np.zeros(npix)
gindices=hp.ang2pix(nside,np.pi/2-B,L)
dgindices=hp.ang2pix(nside,np.pi/2-dB,dL)

ghpmapb=np.zeros(npix)
dghpmapb=np.zeros(npix)
gindicesb=hp.ang2pix(nside,np.pi/2-bB,bL)
dgindicesb=hp.ang2pix(nside,np.pi/2-bdB,bdL)

for sig in range(len(L)):
    ghpmap[gindices[sig]]=ghpmap[gindices[sig]]+1

for back in range(len(bL)):
    ghpmapb[gindicesb[back]]=ghpmapb[gindicesb[back]]+1

for sig in range(len(dL)):
    dghpmap[dgindices[sig]]=dghpmap[dgindices[sig]]+1

for back in range(len(bdL)):
    dghpmapb[dgindicesb[back]]=dghpmapb[dgindicesb[back]]+1

glimamap_l=[] # List for RAW Dataset
dglimamap_l=[] # List for JF12 Deflected Dataset


for k in range(len(dghpmap)): #Use this part for JF12 Deflected Li&Ma projection
    if(dghpmapb[k]==0):
        if(dghpmap[k]==0):
            lima=0
        else:
            print("Error occured : %d"%k)
        dglimamap_l.append(np.nan)
    else:
        lima=LiMa(dghpmap[k],dghpmapb[k]/20,1/20)
        dglimamap_l.append(lima)

hpmap = np.array(dglimamap_l)  
"""
for k in range(len(ghpmap)): #Use this part for RAW Li&Ma projection.
    if(ghpmapb[k]==0):
        if(ghpmap[k]==0):
            lima=0
        else:
            print("Error occured : %d"%k)
        glimamap_l.append(np.nan)
    else:
        lima=LiMa(dghpmap[k],dghpmapb[k]/20,1/20)
        glimamap_l.append(lima)
hpmap = np.array(glimamap_l)
"""   
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

a=np.where(hpmap==np.max(hpmap))
plt.title(r"2D projection of Li&Ma Significance (with Synthetic Injection) : E>5EeV, $N_{side}$=15, $\alpha$=1/20" ,fontsize=28) # title for JF12 backtracked LiMa
plt.show()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
from healpy.newvisufunc import projview

df1=pd.read_csv('Cut2to3EeV.csv')
ID=df1['AugerID'].values
L=df1['L'].values
B=df1['B'].values

df2=pd.read_csv('test2to3EeV.csv')
Trial=df2['Trial'].values
bL=df2['L'].values
bB=df2['B'].values
bdL=df2['dL'].values
bdB=df2['dB'].values

df3=pd.read_csv('JF12_2to3EeVCut.csv')
ID3=df3['AugerID'].values
dL=df3['dL'].values
dB=df3['dB'].values

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

nside=6
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

glimamap_l=[]
dglimamap_l=[]

"""
for k in range(len(dghpmap)):
    if(dghpmapb[k]==0):
        if(dghpmap[k]==0):
            lima=0
        else:
            print("Error occured : %d"%k)
        glimamap_l.append(np.nan)
    else:
        lima=LiMa(dghpmap[k],dghpmapb[k]/20,1/20)
        glimamap_l.append(lima)
"""

for k in range(len(ghpmap)):
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
plt.title("LiMa projection (RAW) : bins=432, alpha=1/20, 2EeV < E $\leq$ 3EeV",fontsize=30)
#plt.figure()
#plt.subplot(111,projection='aitoff')
#plt.hexbin(dl,db)
#plt.hist(hpmap,bins=14)
#plt.yscale('log')
plt.show()

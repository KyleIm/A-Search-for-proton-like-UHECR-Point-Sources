import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp

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
        lima=(n-b)/abs(n-b)*np.sqrt(2*(b/alpha*np.log((b+alpha*b)/(b+alpha*n)))) #Not making inside ln to be 0)
    else:
        lima=(n-b)/abs(n-b)*np.sqrt(2*(n*np.log((n+alpha*n)/(b+alpha*n))+b/alpha*np.log((b+alpha*b)/(b+alpha*n))))

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

glimamap_l=[]
dglimamap_l=[]

print(len(ghpmap),len(dghpmap))
print(dghpmapb[2100]/20, dghpmapb[2150]/20, dghpmapb[2200]/20)
print(dghpmap[2100], dghpmap[2150], dghpmap[2200])

for k in range(len(ghpmap)):
    if(ghpmapb[k]==0):
        if(ghpmap[k]==0):
            lima=0
        else:
            print("g Error occured : %d"%k)
    else:
        lima=LiMa(ghpmap[k],ghpmapb[k]/20,1/20)
        glimamap_l.append(lima)

for k in range(len(dghpmap)):
    if(dghpmapb[k]==0):
        if(dghpmap[k]==0):
            lima=0
        else:
            print("dg Error occured : %d"%k)
            print(dghpmap[k])
    else:
        lima=LiMa(dghpmap[k],dghpmapb[k]/20,1/20)
        dglimamap_l.append(lima)

all_lima_values = np.concatenate((glimamap_l, dglimamap_l))
hist_min = np.min(all_lima_values)
hist_max = np.max(all_lima_values)
if hist_min == hist_max:
    hist_min -= 0.5
    hist_max += 0.5
hist_range = (hist_min, hist_max)


    
plt.figure()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

plt.xlabel("Li-Ma Significance",fontsize = 25)
plt.hist(glimamap_l,bins=30,range=hist_range,facecolor='None',lw=3,edgecolor='blue',label='# of valid bins (RAW) : %d'%len(glimamap_l))
plt.hist(dglimamap_l,bins=30,range=hist_range,facecolor='None',lw=2,edgecolor='red',label='# of valid bins (JF12) : %d'%len(dglimamap_l))
plt.legend(fontsize=20, loc='upper right')
plt.title(r"Histogram of Li&Ma Significance including Synthetic Events,$N_{side}$=15,5EeV,$\alpha$=1/20 - log", fontsize=28)
plt.yscale('log')
plt.show()

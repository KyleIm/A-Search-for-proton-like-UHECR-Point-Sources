import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
from scipy.optimize import curve_fit

df1=pd.read_csv('Cut5EeV.csv')
ID=df1['AugerID'].values
L=df1['L'].values
B=df1['B'].values

df2=pd.read_csv('test5EeV.csv')
Trial=df2['Trial'].values
bL=df2['L'].values
bB=df2['B'].values
bdL=df2['dL'].values
bdB=df2['dB'].values

df3=pd.read_csv('JF12_5EeVCut.csv')
ID3=df3['AugerID'].values
dL=df3['dL'].values
dB=df3['dB'].values

def Gaussian(x,A,mean,std):
    fit_line=A*np.exp(-0.5*((x-mean)/(std*np.sqrt(std)))**2)
    return fit_line

def fit_histogram(values, hist_range, bins=30):
    hist, edges = np.histogram(values, bins=bins, range=hist_range)
    centers = 0.5 * (edges[:-1] + edges[1:])
    mask = hist > 0
    sigma = np.sqrt(hist[mask])
    sigma[sigma == 0] = 1
    p0 = [hist.max(), np.mean(values), np.std(values)]
    params, _ = curve_fit(
        Gaussian,
        centers[mask], hist[mask],
        p0=p0,
        sigma=sigma,
        absolute_sigma=False,
        maxfev=10000
    )
    return params

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

for k in range(len(ghpmap)):
    if(ghpmapb[k]==0):
        if(ghpmap[k]==0):
            lima=0
        else:
            print("Error occured : %d"%k)
    else:
        lima=LiMa(ghpmap[k],ghpmapb[k]/20,1/20)
        glimamap_l.append(lima)

for k in range(len(dghpmap)):
    if(dghpmapb[k]==0):
        if(dghpmap[k]==0):
            lima=0
        else:
            print("Error occured : %d"%k)
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

pp_raw = fit_histogram(glimamap_l, hist_range)
pp_jf = fit_histogram(dglimamap_l, hist_range)


    
plt.figure()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

xfit = np.linspace(hist_min, hist_max, 400)
raw_fit_label = (
    "RAW fit: "
    "$\\mu$={:.4f}, ".format(pp_raw[1]) +
    "$\\sigma$={:.3f}".format(pp_raw[2])
)
jf12_fit_label = (
    "JF12 fit: "
    "$\\mu$={:.4f}, ".format(pp_jf[1]) +
    "$\\sigma$={:.3f}".format(pp_jf[2])
)
plt.plot(xfit, Gaussian(xfit, *pp_raw), lw=3, ls='--', color='blue', label=raw_fit_label)
plt.plot(xfit, Gaussian(xfit, *pp_jf), lw=2, ls='--', color='red', label=jf12_fit_label)

plt.xlabel("Li-Ma Significance",fontsize = 25)
plt.hist(glimamap_l,bins=30,range=hist_range,facecolor='None',lw=3,edgecolor='blue',label='# of valid bins (RAW) : %d'%len(glimamap_l))
plt.hist(dglimamap_l,bins=30,range=hist_range,facecolor='None',lw=2,edgecolor='red',label='# of valid bins (JF12) : %d'%len(dglimamap_l))
plt.legend(fontsize=20, loc='upper right')
plt.title(r"Histogram of Li&Ma Significance: RAW vs JF12,$N_{side}=15$,E>5EeV,$\alpha$=1/20 - log", fontsize=30)
plt.yscale('log')
plt.show()

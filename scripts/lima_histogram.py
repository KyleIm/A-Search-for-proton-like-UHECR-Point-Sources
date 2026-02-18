import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
from scipy.optimize import curve_fit

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

print(len(L),len(dL),len(bL),len(bdL))

def Gaussian(x,A,mean,std):
    fit_line=A*np.exp(-0.5*((x-mean)/(std*np.sqrt(std)))**2)
    return fit_line

def LiMa(n,b,alpha):
    if(b==n):
        lima=np.sqrt(2*(n*np.log((n+alpha*n)/(b+alpha*n))+b/alpha*np.log((b+alpha*b)/(b+alpha*n))))
    elif(n==0):
        lima=(n-b)/abs(n-b)*np.sqrt(2*(b/alpha*np.log((b+alpha*b)/(b+alpha*n)))) #Not making inside ln tobe 0)
    else:
        lima=(n-b)/abs(n-b)*np.sqrt(2*(n*np.log((n+alpha*n)/(b+alpha*n))+b/alpha*np.log((b+alpha*b)/(b+alpha*n))))

    return lima

print("START")

    
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

print("Pass")

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

print("Limaval")
print(len(glimamap_l),len(dglimamap_l))
print(np.min(dglimamap_l),np.min(glimamap_l),np.max(dglimamap_l),np.max(glimamap_l))
print(len(glimamap_l))
H_l,bin_edges=np.histogram(glimamap_l,bins=30,range=(-3.75,3.75))


sigma_l=np.sqrt(H_l)

for j in range(len(sigma_l)):
    if(sigma_l[j]==0):
        sigma_l[j]=1


x_l=np.linspace(-3.75,3.75,30)
ppot,pcov=curve_fit(Gaussian,x_l,H_l,p0=[150,-0.03,0.844],sigma=sigma_l)


gmean=np.mean(glimamap_l)
gstd=np.std(glimamap_l)
gxval=np.linspace(-3.75,3.75,150)
gyval=[]
for i in gxval:
    y=Gaussian(i,ppot[0],ppot[1],ppot[2])
    #y=Gaussian(i,gmean,gstd)
    gyval.append(y)

print(ppot[0],ppot[1],ppot[2])

binsize=int(np.max(dglimamap_l))
H_l,bin_edges=np.histogram(dglimamap_l,bins=30,range=(-3.525,3.755))
sigma_l=np.sqrt(H_l)

for j in range(len(sigma_l)):
    if(sigma_l[j]==0):
        sigma_l[j]=1

dgmean=np.mean(dglimamap_l)
dgstd=np.std(dglimamap_l)
x_l=np.linspace(-3.75,3.75,30)
ppot,pcov=curve_fit(Gaussian,x_l,H_l,p0=[150,-0.03,0.844],sigma=sigma_l)
dgxval=np.linspace(-3.75,3.75,150)
dgyval=[]

print(ppot[0],ppot[1],ppot[2])

for i in dgxval:
    y=Gaussian(i,ppot[0],ppot[1],ppot[2])
    dgyval.append(y)

H_raw, edges = np.histogram(glimamap_l, bins=30, range=(-3.75, 3.75))
centers = 0.5 * (edges[:-1] + edges[1:])
mask = H_raw > 0
p0_raw = [H_raw.max(), np.mean(glimamap_l), np.std(glimamap_l)]
pp_raw, pcov_raw = curve_fit(
    Gaussian,
    centers[mask], H_raw[mask],
    p0=p0_raw,
    sigma=np.sqrt(H_raw[mask]),
    absolute_sigma=False,
    maxfev=10000
)


H_jf, edges_jf = np.histogram(dglimamap_l, bins=30, range=(-3.75, 3.75))
centers_jf = 0.5 * (edges_jf[:-1] + edges_jf[1:])
mask_jf = H_jf > 0
p0_jf = [H_jf.max(), np.mean(dglimamap_l), np.std(dglimamap_l)]
pp_jf, pcov_jf = curve_fit(
    Gaussian,
    centers_jf[mask_jf], H_jf[mask_jf],
    p0=p0_jf,
    sigma=np.sqrt(H_jf[mask_jf]),
    absolute_sigma=False,
    maxfev=10000
)


    
plt.figure()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

xfit = np.linspace(-4, 4, 400)
plt.plot(xfit, Gaussian(xfit, *pp_raw), lw=3, color='blue')
plt.plot(xfit, Gaussian(xfit, *pp_jf), lw=2, color='red')

plt.xlabel("Li-Ma Significance",fontsize = 25)
plt.hist(glimamap_l,bins=30,range=(-4,4),facecolor='None',lw=3,edgecolor='blue',label='RAW : %d'%len(glimamap_l))
plt.hist(dglimamap_l,bins=30,range=(-4,4),facecolor='None',lw=2,edgecolor='red',label='JF12 : %d'%len(dglimamap_l))
plt.legend(fontsize=20)
plt.text(-4.3,ppot[0]*0.3,"A : %0.2f\n$\mu$ : %0.4f \n$\sigma$ : %0.3f"%(ppot[0],ppot[1],ppot[2]),fontsize=25,bbox=dict(boxstyle="square",facecolor='None',edgecolor='black'))
plt.title("Li-Ma : RAW vs JF12(random), bin=432, alpha=1/20, 2EeV < E $\leq$ 3EeV - log", fontsize=30)
plt.yscale('log')
plt.show()

#I need to plot Gaussian Fitting plot

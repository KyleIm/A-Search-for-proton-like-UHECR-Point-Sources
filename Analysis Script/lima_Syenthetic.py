import numpy as np
import pandas as pd
import healpy as hp
from crpropa import *
from scipy.optimize import curve_fit, bisect
from scipy.integrate import quad

df1=pd.read_csv('Cut5EeV.csv')
ID=df1['AugerID'].values
L=df1['L'].values
B=df1['B'].values
E=df1['E'].values

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

Field = JF12Field()
pid = -nucleusId(1, 1)  # charge sign flip for backtracking symmetry
sim = ModuleList()
sim.add(PropagationCK(Field, 1e-8, 0.5 * parsec, 15 * parsec))
obs = Observer()
obs.add(ObserverSurface(Sphere(Vector3d(0), 20 * kpc)))
position = Vector3d(-8.5, 0, 0) * kpc
sim.add(obs)


def LiMa(n,b,alpha):
    if(b==n):
        lima=np.sqrt(2*(n*np.log((n+alpha*n)/(b+alpha*n))+b/alpha*np.log((b+alpha*b)/(b+alpha*n))))
    elif(n==0):
        lima=(n-b)/abs(n-b)*np.sqrt(2*(b/alpha*np.log((b+alpha*b)/(b+alpha*n)))) #Not making inside ln to be 0)
    else:
        lima=(n-b)/abs(n-b)*np.sqrt(2*(n*np.log((n+alpha*n)/(b+alpha*n))+b/alpha*np.log((b+alpha*b)/(b+alpha*n))))

    return lima

def zech_integrand(s, n0, b0):
    n0 = n0
    b0 = b0
    lam = b0 + s
    if lam <= 0:
        return 0.0
    # log-space: lam**n * exp(-lam) / n!  =  exp(n*log(lam) - lam - log(n!))
    log_val = n0 * np.log(lam) - (lam+n0)
    # For extreamly small value put it as zero.
    if log_val < -745:   # np.exp(-745) ~ 5e-324 (lower limit of double)
        return 0.0
    return float(np.exp(log_val))

def zech_upper_limit(n0, b0):
    CL = 0.95
    # Normalize: total integral over [0, ∞)
    Z, _ = quad(zech_integrand, 0, np.inf, args=(n0, b0))

    def integral_diff(sul):
        I, _ = quad(zech_integrand, 0, sul, args=(n0, b0))

        return I / Z - CL

    S_UL = bisect(integral_diff, 0.0, 100.0)
    return S_UL

def E_sample(E):
    return float(np.random.choice(E))

def Deflection(L, B, energy):
    direction = Vector3d()
    lon = L
    lat = np.pi / 2 - B
    direction.setRThetaPhi(1, lat, lon)
    p = ParticleState(pid, energy, position, direction)
    c = Candidate(p)
    sim.run(c)
    d1 = c.current.getDirection()
    dL = d1.getPhi()
    dB = np.pi / 2 - d1.getTheta()
    return dL, dB

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

target_pix = {2100, 2150, 2200}
sul_jf12_map = {}

for i in range(len(ghpmap)):
    if i in target_pix:
        n_raw = ghpmap[i]
        b_raw = ghpmapb[i] / 20
        sul_raw = zech_upper_limit(n_raw, b_raw)
        print("RAW pix %d: n=%d, b=%.6f, zech_upper_limit=%.6f" % (i, n_raw, b_raw, sul_raw))

    if(ghpmapb[i]==0):
        if(ghpmap[i]==0):
            lima=0
        else:
            print("g Error occured : %d"%i)
    else:
        lima=LiMa(ghpmap[i],ghpmapb[i]/20,1/20)
        glimamap_l.append(lima)

for j in range(len(dghpmap)):
    if j in target_pix:
        n_jf12 = dghpmap[j]
        b_jf12 = dghpmapb[j] / 20
        sul_jf12 = zech_upper_limit(n_jf12, b_jf12)
        sul_jf12_map[j] = sul_jf12
        print("JF12 pix %d: n=%d, b=%.6f, zech_upper_limit=%.6f" % (j, n_jf12, b_jf12, sul_jf12))

    if(dghpmapb[j]==0):
        if(dghpmap[j]==0):
            lima=0
        else:
            print("dg Error occured : %d"%j)
    else:
        lima=LiMa(dghpmap[j],dghpmapb[j]/20,1/20)
        dglimamap_l.append(lima)

target_counts = {
    2100: int(np.rint(sul_jf12_map.get(2100, 0.0))),
    2150: int(np.rint(sul_jf12_map.get(2150, 0.0))),
    2200: int(np.rint(sul_jf12_map.get(2200, 0.0))),
}

Sy_2100 = []
Sy_2150 = []
Sy_2200 = []
sy_map = {2100: Sy_2100, 2150: Sy_2150, 2200: Sy_2200}

# 1st-stage rejection set: reject if (L, B, E) already exists in original data.
existing_triplets = set((float(L[k]), float(B[k]), float(E[k])) for k in range(len(E)))
synthetic_triplets = set()

max_trials = 500000
trial = 0

while trial < max_trials:
    trial += 1
    if all(len(sy_map[pix]) >= target_counts[pix] for pix in target_counts):
        break

    # L,B sampled together from the same row; E sampled independently from df1['E'].
    idx = np.random.randint(len(L))
    L_try = float(L[idx])
    B_try = float(B[idx])
    E_try = E_sample(E)

    # Rejection stage 1: reject exact triplet already in original (or already accepted synthetic).
    triplet = (L_try, B_try, E_try)
    if triplet in existing_triplets or triplet in synthetic_triplets:
        continue

    dL_try, dB_try = Deflection(L_try, B_try, E_try)
    pix_try = hp.ang2pix(nside, np.pi/2 - dB_try, dL_try)

    # Rejection stage 2: reject if not in target indices or target already filled.
    if pix_try not in target_counts:
        continue
    if len(sy_map[pix_try]) >= target_counts[pix_try]:
        continue

    sy_map[pix_try].append((L_try, B_try, dL_try, dB_try, E_try))
    synthetic_triplets.add(triplet)

print("target_counts =", target_counts)
print("filled_counts =", {2100: len(Sy_2100), 2150: len(Sy_2150), 2200: len(Sy_2200)})
print("trials =", trial)

output_path = "Sy_events_5EeV.csv"
with open(output_path, "w") as f:
    f.write("Indices,L,B,dL,dB,E\n")
    for pix in [2100, 2150, 2200]:
        for L_val, B_val, dL_val, dB_val, E_val in sy_map[pix]:
            line = "%d, %lf, %lf, %lf, %lf, %lf\n" % (pix, L_val, B_val, dL_val, dB_val, E_val)
            f.write(line)

print("saved_csv =", output_path)

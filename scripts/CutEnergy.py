import numpy as np

rf=open('Archive_v6r2p2', 'r')
wf=open('Cut50EeV.csv','w')
head="AugerID,Theta,Phi,RA,Dec,L,B,E,UTC,GPStime\n"
wf.write(head)

E1=0
E2=0
E3=0
E4=0
E5=0

n=0
for i in rf:
    par=i.split(' ')
    ID=int(par[0])
    if(len(par)>36):
        E=float(par[37])
        Theta=np.deg2rad(float(par[2]))
        Phi=np.deg2rad(float(par[3]))
        L=np.deg2rad(float(par[5]))
        B=np.deg2rad(float(par[6]))
        UTC=int(par[7])
        RA=np.deg2rad(float(par[13]))
        Dec=np.deg2rad(float(par[14]))
        GPST=int(par[39])
        par22=int(par[21])
        par23=int(par[22])
        par43=int(par[42])
        par44=int(par[43])
        if((par22>0 and par23>1 and Theta<np.pi/3) and E>50 and par43>3 and (par44==1 or par44==2)): #Required conditions for cut
            data="%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d\n"%(ID,Theta,Phi,RA,Dec,L,B,E,UTC,GPST)
            wf.write(data)
            n=n+1
print(n)
#print(E1,E2,E3,E4,E5)
rf.close()
wf.close()

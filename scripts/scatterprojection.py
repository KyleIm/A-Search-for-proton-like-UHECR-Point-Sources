import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
from healpy.newvisufunc import projview

df=pd.read_csv('JF12_50EeVCut.csv')
#ID=df['AugerID'].values
#E=df['E'].values
l=df['L'].values * -1
b=df['B'].values


plt.figure()
ax2=plt.subplot(111, projection="mollweide")
plt.scatter(l,b,c='b',s=10)
ax=plt.gca()
ax.set_xticklabels(ax.get_xticklabels()[::-1])
ax2.tick_params(axis='both', which='major', labelsize=20)
plt.title("Selected Event location for magnetic field", fontsize = 30)
plt.xlabel("L", fontsize = 25)
plt.ylabel("B", fontsize = 25)
plt.grid(True)

plt.show()

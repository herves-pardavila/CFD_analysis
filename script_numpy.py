# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 14:25:56 2025

@author: e6040
"""

import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
plt.close("all")

main_path='C:/TFM_DAVIDH_2025/'

datos=np.loadtxt(main_path+"Tablas/"+"ustar0.33_NO_7200_PS_T_V_threshold_plane_z=3m.csv",delimiter=",",skiprows=1)

x=datos[:,-3]
y=datos[:,-2]
area=datos[:,4]

ps1=datos[:,0]
ps2=datos[:,1]
ps3=datos[:,2]


# #interpolación, tarda mucho con 8GB de RAM y no aporta mucho
# minx,maxx,dimx=min(x),max(x),len(x)
# miny,maxy,dimy=min(y),max(y),len(y)

# xt=np.linspace(minx,maxx,dimx)
# yt=np.linspace(miny,maxy,dimy)

# X,Y=np.meshgrid(xt,yt)


# Z = griddata((x, y), ps1, (X,Y), method='nearest')
# Z[ ( (X > 115) & (X < 215) ) & (Y > 115) & (Y < 215)]=np.nan

# fig=plt.figure()
# plt.pcolormesh(X,Y,Z,cmap="coolwarm")
# plt.savefig("imagen_gridded_numpy.png")
# plt.show()

#scatter plot

fig=plt.figure()
ax=fig.add_subplot(111)
im=ax.scatter(x,y,s=area**2,c=ps1,marker="h",cmap="coolwarm",edgecolor=None)
ax.set_aspect('equal', adjustable='box')
ax.set_facecolor('black')
fig.colorbar(im,location="top")
fig.savefig("imagen_scatter_numpy.png")
plt.show()

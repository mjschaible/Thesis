
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate

#while cont != 0:
#    n = n+1
#filename = raw_input('Enter filename: ')
filename1 = 'time000all.d'

contents = np.genfromtxt(filename1, unpack=True)

xcoord = contents[0,:]
ycoord = contents[1,:]
zcoord = contents[2,:]
    
temp = contents[3,:]
dens = contents[4,:]
conc = contents[5,:]

'''
c_temp = np.array(temp)
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111, projection='3d')
cmhot = plt.get_cmap("rainbow")
cax = ax1.scatter(xcoord, ycoord, zcoord, s = 10, c = c_temp, cmap=cmhot, lw=0)
'''

filename2 = 'Contour.out'
contents = np.genfromtxt(filename2, unpack=True)

time = contents[0,:]
zpos = contents[1,:]
avgtemp = contents[2,:]
    
dens = contents[3,:]
conc = contents[4,:]

cols = np.unique(zpos).shape[0]
Time = time.reshape(-1, cols)
Zpos = zpos.reshape(-1, cols)
ATemp = avgtemp.reshape(-1,cols)
#print cols, Time

fig3 = plt.figure(3)
ax3 = fig3.add_subplot(111)
ax3.imshow(ATemp.T, origin='lower', extent=[time.min(), time.max(), zpos.min(), zpos.max()])
ax3.set_aspect(0.25)
#plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
#plt.pcolormesh(xi,yi,zi,cmap=cmhot)

'''
fig2=plt.figure(2)
ax2=fig2.add_subplot(111)
c_temp=np.array(avgtemp)
cax = ax2.scatter(time,zpos,c=c_temp,cmap=cmhot, lw=0)
#ax2.colorbar(shrink=0.92)
'''
plt.show()

#import numpy as np
#import scipy 
from pylab import *
from scipy.optimize import curve_fit

cem_vf = np.array([1e-7, 1e-6, 1e-5, 1e-4])
Rcon = np.array([0.66, 1.187, 2.109, 3.752])

kcem65= np.array([0.161, 0.283, 0.48, 0.808])
kcem20= np.array([0.155, 0.268, 0.443, 0.721])
kcem10= np.array([0.152, 0.253, 0.408, 0.646])
kcem01= np.array([0.107, 0.160, 0.223, 0.299])

def fitexp(x, a, c, d):
    return a*np.exp(-c*x)+d

def fitpower(x, m, q):
    return q*(x**m)


#popt65, pcov20 = curve_fit(fitexp, Rcon/25, kcem65, p0 = (5, 1, 1))
#popt20, pcov65 = curve_fit(fitexp, Rcon/25, kcem20, p0 = (5, 1, 1))

popt65, pcov65 = curve_fit(fitpower, Rcon/25, kcem65, p0 = (2, 5))
popt20, pcov20 = curve_fit(fitpower, Rcon/25, kcem20, p0 = (2, 5))
popt10, pcov10 = curve_fit(fitpower, Rcon/25, kcem10, p0 = (2, 5))
popt01, pcov01 = curve_fit(fitpower, Rcon/25, kcem01, p0 = (2, 5))
print popt65[0]
print popt20[0]
print popt10[0]
print popt01[0]

#xx65=np.logspace(-8,-2,num=10)
#yy65 = fitexp(xx65, *popt65)
#xx20=np.logspace(-8,-2,num=10)
#yy20 = fitexp(xx20, *popt20)
#print xx, yy
xx65=np.logspace(-2,1,num=1000)
yy65 = fitpower(xx65, *popt65)
xx20=np.logspace(-2,1,num=1000)
yy20 = fitpower(xx20, *popt20)
xx10=np.logspace(-2,1,num=1000)
yy10 = fitpower(xx10, *popt10)
xx01=np.logspace(-2,1,num=1000)
yy01 = fitpower(xx01, *popt01)


print xx65[0], yy65[0]
check = 4.59771172*0.01**(0.91602145)
print check

fig = plt.figure(1)
ax1 = fig.add_subplot(111)
ax1.semilogx(Rcon/25, kcem65, 'ko')
ax1.semilogx(Rcon/25, kcem20, 'ro')
ax1.semilogx(Rcon/25, kcem10, 'bo')
ax1.semilogx(Rcon/25, kcem01, 'go')
ax1.semilogx(xx65, yy65, linestyle = '-', label = 'kcem=6.5', color='k')
ax1.semilogx(xx20, yy20, linestyle = '-', label = 'kcem=2.0', color='r')
ax1.semilogx(xx10, yy10, linestyle = '-', label = 'kcem=1.0', color='b')
ax1.semilogx(xx01, yy01, linestyle = '-', label = 'kcem=0.1', color='g')
ax1.set_xlim(1e-2, .2)
ax1.set_ylim(0,1)
ax1.set_xlabel('Contact to Grain Radius Ratio (Rcon/rg)')
ax1.set_ylabel('Effective Thermal Conductivity [W/mK]')

ax1.text(.011, .75, 'Fit Parameters:', fontsize=10)
ax1.text(.011, .7, 'kcem=6.5: c*x^{%1.2f}' % popt65[0], fontsize=10)
ax1.text(.011, .65, 'kcem=2.0: c*x^{%1.2f}' % popt20[0], fontsize=10)
ax1.text(.011, .6, 'kcem=1.0: c*x^{%1.2f}' % popt10[0], fontsize=10)
ax1.text(.011, .55, 'kcem=0.1: c*x^{%1.2f}' % popt01[0], fontsize=10)
ax1.legend(loc = 'upper left', prop={'size':10})

print 'kcem=0.1: c*x^{0}'.format(popt01[0])
#ax1.set_title('Wood Linear')


plt.savefig('./piqueuxFit',dpi=600)

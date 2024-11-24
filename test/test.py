import numpy as np
import okada4py as ok92

#-----------------------------------
# 1. Make a list of receivers
#-----------------------------------

xs = np.arange(-20.0, 20.0, 0.1)
ys = np.arange(-20.0, 20.0, 0.1)
xs, ys = np.meshgrid(xs, ys)
xs = xs.flatten()
ys = ys.flatten()
zs = np.zeros(xs.shape)
zrec = zs + 5*(np.exp(-1*((xs/4+2)**2 + (ys/4+2)**2))) #topographic surface

#-----------------------------------
# 2. Make a dislocation
#-----------------------------------

xc = np.array([0.0])
yc = np.array([0.0])
depth = np.array([2.12133])
length = np.array([10.0])
width = np.array([6.0])
dip = np.array([45.])
strike = np.array([0.0])

ss = np.array([1.0])
ds = np.array([0.0])
ts = np.array([0.0])

mu = 30.0e9
nu = 0.25

#-----------------------------------
# 3. Run the routine for okada 92
#-----------------------------------

u, d, s, flag, flag2 = ok92.okada92(xs, ys, zs, xc, yc, depth, length, width, dip, strike, ss, ds, ts, mu, nu)
u = u.reshape((xs.shape[0], 3))
s = s.reshape((xs.shape[0], 6))
stress = np.zeros((3, 3, len(xs)))
stress[0,0,:] = s[:,0]
stress[1,1,:] = s[:,3]
stress[2,2,:] = s[:,5]
stress[0,1,:] = s[:,1]
stress[1,0,:] = s[:,1]
stress[0,2,:] = s[:,2]
stress[2,0,:] = s[:,2]
stress[1,2,:] = s[:,4]
stress[2,1,:] = s[:,4]
Trace = np.trace(stress, axis1=0, axis2=1)
a = np.tile(np.eye(3), (len(xs), 1, 1)).T
a[0,0,:] = Trace
a[1,1,:] = Trace
a[2,2,:] = Trace
stress = stress - a
Frob = np.array([np.linalg.norm(stress[:,:,i]) for i in range(len(xs))])
Sigmaxy = stress[0,1,:]

#-----------------------------------
# Finite difference
#-----------------------------------

x = xs.reshape((400,400))
y = ys.reshape((400,400))

Ux = u[:,0]
Ux = Ux.reshape((400,400))

Uy = u[:,1]
Uy = Uy.reshape((400,400))

gradientx = np.gradient(Ux, 0.1)
gradienty = np.gradient(Uy, 0.1)

Uxx = gradientx[1]
Uxy = gradientx[0]
Uyx = gradienty[1]
Uyy = gradienty[0]

Uxy = 0.5*(Uxy + Uyx)
Sxy = mu*Uxy

#------------------------------------------------
# 4. Run the okada92 routine with topographic DEM
#------------------------------------------------
ut, dt, st, flagt, flag2t = ok92.okada92(xs, ys, zs, xc, yc, depth, length, width, dip, strike, ss, ds, ts, mu, nu, zrec)
ut = ut.reshape((xs.shape[0], 3))


#-----------------------------------
# 5. Plot 
#-----------------------------------

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('image', cmap='seismic')

plt.figure(1)
plt.suptitle('Okada 92')
plt.subplot(131)
plt.scatter(xs, ys, s=1.0, c=u[:,0], linewidth=0.0) 
plt.colorbar(orientation='horizontal', shrink=1.)
plt.subplot(132)
plt.scatter(xs, ys, s=1.0, c=u[:,1], linewidth=0.0) 
plt.colorbar(orientation='horizontal', shrink=1.)
plt.subplot(133)
plt.scatter(xs, ys, s=1.0, c=u[:,2], linewidth=0.0) 
plt.colorbar(orientation='horizontal', shrink=1.)

plt.figure(2)
plt.title('Stress')
plt.scatter(xs, ys, s=1.0, c=Sigmaxy, linewidth=0.0, vmin=-1e9, vmax=1e9)
plt.colorbar()

plt.figure(3)
plt.title("Topographic Surface")
plt.scatter(xs, ys, c=zrec, linewidth=0.0, cmap='terrain')
plt.colorbar()

plt.figure(4)
plt.suptitle('Okada 92 with topographic correction')
plt.subplot(131)
plt.scatter(xs, ys, s=1.0, c=ut[:,0], linewidth=0.0) 
plt.colorbar(orientation='horizontal', shrink=1.)
plt.subplot(132)
plt.scatter(xs, ys, s=1.0, c=ut[:,1], linewidth=0.0) 
plt.colorbar(orientation='horizontal', shrink=1.)
plt.subplot(133)
plt.scatter(xs, ys, s=1.0, c=ut[:,2], linewidth=0.0) 
plt.colorbar(orientation='horizontal', shrink=1.)

plt.show()

from boutdata import collect
from boututils.showdata import showdata
import matplotlib.pyplot as plt
import numpy as np

phi = collect("total_phi")
phiext = collect("phiext")

nt,nx,ny,nz = phi.shape

ne = collect("Ne")
ni = collect("Ni")
ng = collect("Ng")
NeE = collect("NeE")
n0 = collect("n0")
ne = ne[:,2:-2,:,:]
ni = ni[:,2:-2,:,:]
ng = ng[:,2:-2,:,:]
NeE = NeE[:,2:-2,:,:]
phi = phi[:,2:-2,:,:]
ne *= n0
ni *= n0
ng *= n0
NeE *= n0

nt,nx,ny,nz = phi.shape
print(phi.shape)

t = collect("t_array")
w0 = collect("w0")
t /= w0 / 1e6

dx = collect("dx")
L0 = collect("L0")
dx = dx[0,0]
xx = np.linspace(0,dx*nx,nx)
xx *= L0

print(xx.shape)
print(ni.shape)

# plot densities
plt.figure(figsize=(15,5))
plt.subplot(131)
plt.title('electron density')
plt.contourf(xx,t,ne[:,:,0,0],100)
plt.xlabel('position')
plt.ylabel(r'time [$\mu$s]')
plt.colorbar()
plt.subplot(132)
plt.title('ion density')
plt.contourf(xx,t,ni[:,:,0,0],100)
plt.xlabel('position')
plt.ylabel(r'time [$\mu$s]')
plt.colorbar()
plt.subplot(133)
plt.title('neutral gas density')
plt.contourf(xx,t,ng[:,:,0,0],100)
plt.xlabel('position')
plt.ylabel(r'time [$\mu$s]')
plt.colorbar()
plt.tight_layout()
plt.show()

# plot energy density
plt.figure(figsize=(5,5))
plt.subplot(111)
plt.title('electron energy density')
plt.contourf(xx,t,NeE[:,:,0,0],100)
plt.xlabel('position')
plt.ylabel(r'time [$\mu$s]')
plt.colorbar()
plt.tight_layout()
plt.show()

showdata(ne[::2,:,0,:])
showdata(ni[::2,:,0,:])
showdata(ng[::2,:,0,:])

showdata(phi[:,:,0,:])

plt.plot(phiext[:,nx//2,0,0])
plt.grid()
plt.show()

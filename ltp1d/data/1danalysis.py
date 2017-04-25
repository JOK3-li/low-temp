from boutdata import collect
from boututils.showdata import showdata
import matplotlib.pyplot as plt
import numpy as np

phi = collect("total_phi")
phiext = collect("phiext")

nt,nx,ny,nz = phi.shape

tstart = 0

ne = collect("Ne")
ni = collect("Ni")
# ng = collect("Ng")
NeE = collect("NeE")
n0 = collect("n0")
T0 = collect("T0")
phi0 = collect("v0")
ne = ne[tstart:,2:-2,:,0]
ni = ni[tstart:,2:-2,:,0]
# ng = ng[:,2:-2,:,0]
NeE = NeE[tstart:,2:-2,:,0]
phi = phi[tstart:,2:-2,:,0]
E = NeE / ne
ne *= n0
ni *= n0
# ng *= n0
NeE *= n0*T0
E *= T0 # eV
phi *= phi0
phiext *= phi0

nt,nx,ny = phi.shape

t = collect("t_array")
w0 = collect("w0")
t /= w0 / 1e6
t = t[tstart:]

dx = collect("dx")
print(dx.shape)
L0 = collect("L0")
dx *= L0
xx = np.cumsum(dx[2:-2,0])

# plot densities
plt.figure(figsize=(15,5))
plt.subplot(131)
plt.title('electron density')
plt.contourf(xx,t,ne[:,:,0],100)
plt.xlabel('position')
plt.ylabel(r'time [$\mu$s]')
plt.colorbar()
plt.subplot(132)
plt.title('ion density')
plt.contourf(xx,t,ni[:,:,0],100)
plt.xlabel('position')
plt.ylabel(r'time [$\mu$s]')
plt.colorbar()
plt.subplot(133)
plt.title('electron energy density')
plt.contourf(xx,t,np.log(E[:,:,0]),100)
plt.xlabel('position')
plt.ylabel(r'time [$\mu$s]')
plt.colorbar()
plt.tight_layout()
plt.show()

showdata(ne[::2,:,0])
showdata(ni[::2,:,0])
# showdata(ng[::2,:,0])
showdata(E[:,:,0])
showdata(phi[:,:,0])

plt.plot(phiext[:,-1,0])
plt.grid()
plt.show()

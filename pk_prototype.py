#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

"""
pk_prototype: code to calculate Power Spectrum from a N-body simulation's (with Periodic BC)
particles or halos.

"""


"""
FIRST: get count in cells grid

Simple version: use only x,y,z, no weights yet

Steps:
- Read the file 
- Define grid size (try Nx=Ny=Nz< 128)
- process = assign densities to grid
- output: an array w/densities.. to be used by other program....

work:
1) Read a file and plot points (with a fixed file as of now) 
2) Make the file a variable

today is Wed 3/4/2013
"""

import numpy as np
import matplotlib.pyplot as plt


#load figure (  is there an smater way to do this?)
a = np.loadtxt('ghalos_z0p0_dlow_L8_01_real.dat.xyz') #change later to a variable file
x = np.array(a[:,0])
y = np.array(a[:,1])
z = np.array(a[:,2])

print "loaded file... \n"

#plot points
plt.figure()
plt.plot(x,y,'.')
plt.grid()
plt.show()
print "plot points...\n"

#get cell...
Lbox = 500.0
Lx = Ly = Lz = Lbox
Vbox = Lx*Ly*Lz
Nx = Ny = Nz = 64
rho = np.zeros((Nx,Ny,Nz))

n_x = (np.floor((x/Lx)*Nx))
n_y = (np.floor((y/Ly)*Ny))
n_z = (np.floor((z/Lz)*Nz))

print "calculated cell postions...\n"

for i in range(len(n_x)): #cloud in cell stuff
    dx = (Nx/Lx)*x[i] - n_x[i] + 0.5
    tx = 1 - abs(dx)
    dy = (Ny/Ly)*y[i] - n_y[i] + 0.5
    ty = 1 - abs(dy)
    dz = (Nz/Lz)*z[i] - n_z[i] + 0.5
    tz = 1 - abs(dz)
    #cloud in cell densities
    i1 = (int)(n_x[i])
    i2 = (int)(np.mod(i1 + np.sign(dx),Nx))
    j1 = (int)(n_y[i])
    j2 = (int)(np.mod(j1 + np.sign(dy),Ny))
    k1 = (int)(n_z[i])
    k2 = (int)(np.mod(k1 + np.sign(dy),Nz))

    rho[i1,j1,k1] += tx*ty*tz/(Vbox/Nx/Ny/Nz)
    rho[i2,j1,k1] += dx*ty*tz/(Vbox/Nx/Ny/Nz)
    rho[i1,j2,k1] += tx*dy*tz/(Vbox/Nx/Ny/Nz)
    rho[i2,j2,k1] += dx*dy*tz/(Vbox/Nx/Ny/Nz)
    rho[i1,j1,k2] += tx*ty*dz/(Vbox/Nx/Ny/Nz)
    rho[i2,j1,k2] += dx*ty*dz/(Vbox/Nx/Ny/Nz)
    rho[i1,j2,k2] += tx*dy*dz/(Vbox/Nx/Ny/Nz)
    rho[i2,j2,k2] += dx*dy*dz/(Vbox/Nx/Ny/Nz)

#transform to overdensities

rhomean = len(n_x)/(Lbox**3.0)
print "my mean= "+str(rhomean)
print "true mean= "+str(rho.mean())
delta = rho/rhomean - 1.0;
print "delta mean= "+str(delta.mean())
a_forplot = np.zeros((Nx,Ny))


for i in range(10):
    a_forplot +=  delta[:,:,30+i]

plt.figure()
plt.imshow(a_forplot, origin = 'lower')
plt.grid()
plt.show()
print delta.min()
print delta.max()

#now, get the fourier transform of delta

#in scipy
#import scipy

delta_k = np.fft.fftn(delta)

#now get the raw power spectrum

#simplest way for now:
# for all i,j,half of k: find out if total k falls into bin, count it in

tot_fk = 0.0;
tot_kmodes = 0.0;

for ix in range(Nx):
    if (ix<=Nx/2):
        kx = 2.0*np.pi*ix/Lx
    else:
        kx = 2.0*np.pi*(ix - Nx)/Lx
    for jy in range(64):
        if (jy<=Ny/2):
            ky = 2.0*np.pi*jy/Ly   
        else:
            ky = 2.0*np.pi*(jy - Ny)/Ly
        for lz in range(33):
            kz = 2.0*np.pi*lz/Lz
            ktot = np.sqrt(kx*kx + ky*ky + kz*kz)
            if (ktot > 0.1 and ktot < 0.2):
                tot_fk = tot_fk +  delta_k[ix,jy,lz]*delta_k[-ix,-jy,-lz] #rectify this later
                tot_kmodes = tot_kmodes + 1

print tot_fk
print tot_kmodes












                


                





    















#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python


""" cic_simple.py: Module to obtain the cloud-in-cell density grid from a set of points
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

#l = x[0]
#print l
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

#print len(n_x)
#print n_x[50]
#print n_y[50]
#print n_z[50]

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


"""

#now, calculate two-point function
print "calculate 2PCF...\n"
print "find compatible cells...\n"
r_xi = 30.0 #Mpc/h
dr_xi = 5.0 #don't worry about \Delta_r for now
r_ximin = r_xi - dr_xi/2.0
r_ximax = r_xi + dr_xi/2.0
# a particle in the middle of the mesh (at Nx/2,Ny/2,Nz/2) has - VERIFY
x_cell = Lbox/(2*Nx) + 0.5*Lbox;
y_cell = x_cell;
z_cell = x_cell;
x_yes = np.zeros(0)
y_yes = np.zeros(0)
z_yes = np.zeros(0)

cell_yes = 0
for i in range(Nx):
    x_pos = Lbox/(2*Nx) + i*Lbox/(Nx+0.0)
    diff_x = (x_cell - x_pos)**2
    for j in range(Ny):
        y_pos = Lbox/(2*Ny) + j*Lbox/(Ny+0.0)
        diff_y = (y_cell - y_pos)**2
        for k in range(Nz/2):#change to be on top semi-sphere
            z_pos = Lbox/(2*Nz) + k*Lbox/(Ny+0.0)
            diff_z = (z_cell - z_pos)**2
            dist_cells = np.sqrt(diff_x + diff_y + diff_z)
            if (dist_cells > r_ximin):
                if(dist_cells < r_ximax):
                    x_yes=np.append(x_yes, (i - Nx/2))
                    y_yes=np.append(y_yes, (j - Ny/2))
                    z_yes=np.append(z_yes, (k - Nz/2))
                    cell_yes = cell_yes + 1



x_yes = np.int_(x_yes)
y_yes = np.int_(y_yes)
z_yes = np.int_(z_yes)
#now get the pairs
print "now get the pairs...\n"
dd = 0.0
pairs = 0

for i in range(Nx):
    print i
    id_2 = np.mod(i+x_yes,Nx)
   # print id_2
   # print np.amin(id_2)
   # print np.amax(id_2)
    for j in range(Ny):
        jd_2 = np.mod(j+y_yes,Ny)
        for k in range(Nz):#change to be on top semi-sphere
            delta1 = delta[i][j][k]
            kd_2 = np.mod(k+z_yes,Nz)
            #np.min(kd_2)
            #np.max(kd_2)
            for c in range(cell_yes):
                delta2 = delta[id_2[c]][jd_2[c]][kd_2[c]]
                dd = dd + delta1*delta2
                pairs = pairs + 1

        
print dd
print pairs
print dd/pairs

print "done.\n"

"""







                


                





    















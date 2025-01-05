import sys #https://scipython.com/blog/visualizing-the-earths-magnetic-field/
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# Mean magnitude of the Earth's magnetic field at the equator in T
B0 = 3.12e-5
# Radius of Earth, Mm (10^6 m: mega-metres!)
RE = 6.370
# Deviation of magnetic pole from axis
alpha = np.radians(9.6)

def B(r, theta):
    """Return the magnetic field vector at (r, theta)."""
    fac = B0 * (RE / r)**3
    return -2 * fac * np.cos(theta + alpha), -fac * np.sin(theta + alpha)

# Grid of x, y points on a Cartesian grid
nx, ny, nz = 500, 500, 500
XMAX, YMAX, ZMAX = 40, 40, 40
x = np.linspace(-XMAX, XMAX, nx)
y = np.linspace(-YMAX, YMAX, ny)
z = np.linspace(-ZMAX, ZMAX, nx)
X, Y = np.meshgrid(x, y)
r, theta = np.hypot(X, Y), np.arctan2(Y, X)

# Magnetic field vector, B = (Ex, Ey), as separate components
Br, Btheta = B(r, theta)
# Transform to Cartesian coordinates: NB make North point up, not to the right.
c, s = np.cos(np.pi/2 + theta), np.sin(np.pi/2 + theta)
Bx = -Btheta * s + Br * c
By = Btheta * c + Br * s
Bz= Bx
Bxyz=np.array([Bx,By,Bz])
fig, ax = plt.subplots()

# Plot the streamlines with an appropriate colormap and arrow style
color = 2 * np.log(np.hypot(Bx, By))
ax.streamplot(x, y, Bx, By, color=color, linewidth=1, cmap=plt.cm.inferno, density=0.4, arrowstyle='->', arrowsize=1.5)

# Add a filled circle for the Earth; make sure it's on top of the streamlines.
ax.add_patch(Circle((0,0), RE, color='b', zorder=100))

def motion():
    m = 9.11E-31
    q = -1.6E-19
    QoverM = q/m

    dt = 0.1 #small timestep

    t = np.arange(0.0, 100.0, dt) #create an array that will hold the times
    rp = np.zeros((len(t), 3)) #create an array that will hold the position values
    vp = np.zeros((len(t), 3)) #create an array that will hold the velocity values

    v0 = 100.0 #set the initial velocity to 100 m/s
    rp[0,:] = np.array([20, 0, 0]) #initialize the position to y=-5, 5m above the lower dipole
    vp[0,:] = np.array([100, 0, 0]) #initialize the velocity to be in the z-direction
    print(rp)
    
    for it in np.arange(0, len(t)-1,1):
        rq, thetaq = np.hypot(rp[it,0], rp[it, 1]), np.arctan2(rp[it,1], rp[it, 0])

        # Magnetic field vector, B = (Ex, Ey), as separate components
        Brq, Bthetaq = B(rq, thetaq)
        c, s = np.cos(np.pi/2 + thetaq), np.sin(np.pi/2 + thetaq)
        Bp=[ -Bthetaq * s + Brq * c, Bthetaq * c + Brq * s, -Bthetaq * s + Brq * c]
        
        Ap = QoverM * np.cross(vp[it,:], Bp) #Calculate the magnetic force on the particle
        vp[it+1] = vp[it] + dt*Ap #Update the velocity of the particle based on this force
        rp[it+1] = rp[it] + dt*vp[it] #Update the positon of the particle based on this velocity
        #if (np.sqrt(np.sum(rp[it+1]**2)) > 400.0): #If the particle escapes (goes more than 20m away from the origin) end the loop
           # break
    return rp
rp = motion()
print(rp)
# now to make different views of the charged particle's trajectory
#plt.streamplot(Y,Z, Bxyz[:,:,1], Bxyz[:,:,2], color="black")
plt.plot(rp[:,0], rp[:,1])


ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_xlim(-XMAX, XMAX)
ax.set_ylim(-YMAX, YMAX)
ax.set_aspect('equal')
plt.show()
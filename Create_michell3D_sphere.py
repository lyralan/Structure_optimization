import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from numpy import linalg as LA
from mpl_toolkits.mplot3d import axes3d, Axes3D
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--r', type=int, default=30)
parser.add_argument('--r1', type=int, default=5)
parser.add_argument('--r2', type=int, default=10)
parser.add_argument('--alpha', type=float, default=8)
args = parser.parse_args()

def sph2cart(azimuth,elevation,r):
    x = r * np.cos(elevation) * np.cos(azimuth)
    y = r * np.cos(elevation) * np.sin(azimuth)
    z = r * np.sin(elevation)
    return x, y, z

def michell_3d_sphere(r,r1,r2,alpha): 
    # step 1 Domain     
    step = 2                                        # step in absolute value
    # step 2 Loads
    m = np.array([[0,0,0],
                  [0,0,1],
                  [0,1,0]])                         # strain tensor

    # step 3 Generation of the form
    c1 = np.array([np.sqrt(r**2-r1**2), 0, 0])      # Position of Torque 1 center
    c2 = np.array([np.sqrt(r**2-r2**2), 0, np.pi])  # Position of Torque 2 center

    # step 4 Define starting points on Torque 1
    starting_points = [] #  n x 3
    long = 2*np.pi/alpha

    for i in 1 + np.arange(long):
        a = 2 * np.pi * i / long
        node = np.array([r,a,np.arccos(c1[0]/r)])
        starting_points.append(node)
    starting_points = np.array(starting_points)

    #step 5 Calculate eigenvectors and eigenvalues
    d, v = LA.eig(m)
    v = v[:, d != 0]
    d = d[d != 0]
    lambda_M = d[0]
    lambda_m = d[-1]
    pi_M = v[:, 0] * np.sign(v[-1, 0])
    pi_m = v[:, -1] * np.sign(v[-1, -1])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    next_starting_points = []

    # step 6 First iteration
    for r, a, angle in starting_points:
        x, y, z = sph2cart(a, np.pi/2-angle, r)
        for (pi_r, pi_a, pi_angle), color in (pi_M, 'dodgerblue'), (pi_m, 'hotpink'):
            mr, ma, mangle = m_sphere = [r + step * pi_r, 
                                         a + step * pi_a/(r * np.sin(angle)), 
                                         angle + step * pi_angle/r] 
            mx, my, mz = sph2cart(ma, np.pi/2 - mangle, mr)
            ax.plot([x, mx],
                    [y, my],
                    zs = [z, mz], 
                    color = color,
                    linewidth = 0.5)
            next_starting_points.append(m_sphere)

    # step 7 Next Iterations
    while max(np.array(next_starting_points)[:,2]) < (np.pi-np.arccos(c2[0]/r)):
            starting_points = next_starting_points
            next_starting_points = []
            
            for r, a, angle in starting_points:
                x, y, z = sph2cart(a, np.pi/2 - angle, r)
                time = len(next_starting_points) + 1
                pi_select, color = (pi_M, 'dodgerblue') if time % 2 == 1 else (pi_m, 'hotpink')    
                [pi_r, pi_a, pi_angle] = pi_select
                mr, ma, mangle = m_sphere = [r + step * pi_r, 
                                             a + step * pi_a/(r * np.sin(angle)), 
                                             angle + step * pi_angle/r] 
                mx, my, mz = sph2cart(ma, np.pi/2 - mangle, mr)
                ax.plot([x, mx],
                        [y, my],
                        zs = [z, mz], 
                        color = color,
                        linewidth = 0.5)
                next_starting_points.append(m_sphere)
    for ii in np.arange(0,360,30):
        ax.view_init(elev=90., azim=ii)
        fig.savefig("movie%d.png" % ii)
    plt.show()
michell_3d_sphere(args.r, args.r1, args.r2, np.pi / args.alpha) #user input radius of sphere, torque 1, torque 2, radial distribution of starting points 

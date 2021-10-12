import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from numpy import linalg as LA
from mpl_toolkits.mplot3d import axes3d, Axes3D
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--point1', type=int, default=16)
parser.add_argument('--point2', type=int, default=4)
args = parser.parse_args()

def sph2cart(azimuth,elevation,r):
        x = r * np.cos(elevation) * np.cos(azimuth)
        y = r * np.cos(elevation) * np.sin(azimuth)
        z = r * np.sin(elevation)
        return x, y, z

def michell_3d_spinning(point1, point2):       
    # Domain
    r1, r2, step = 10, 35, 1        # Radius of Bowls 1 & 2, Step in absolute value
    hstep = 0                       # For exploded view

    # Loads
    m = np.array([[0,1,0],
                  [1,0,0],
                  [0,0,0]])         # Strain tensor for this coordinate sys. convention

    # Calculate eigenvectors and eigenvalues
    d, v = LA.eig(m)
    v = v[:, d != 0]
    d = d[d != 0]
    lambda_M = d[0]
    lambda_m = d[-1]
    pi_M = v[:, 0] * np.sign(v[0, 0])
    pi_m = v[:, -1] * np.sign(v[0, -1])

    # Generation of the Form
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = r1 * np.outer(np.cos(u), np.sin(v))
    y = r1 * np.outer(np.sin(u), np.sin(v))
    z = r1 * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color=(0.8, 0.8, 0.8), alpha=0.4)

    for i in np.arange(1,point2):
        phi = np.pi * i / point2
        starting_points = []
        for ii in 1 + np.arange(point1 * 2):
            th = 2 * np.pi * ii / (point1 * 2)   
            node = np.array([r1, th, phi])
            starting_points.append(node)
        next_starting_points = []    
        height = (1 - phi / (np.pi / 2)) * hstep 
        
        # First iteration
        for r, a, angle in starting_points:
            x, y, z = sph2cart(a, np.pi/2 - angle, r)
            for (pi_r, pi_a, pi_angle), color in (pi_M, 'dodgerblue'), (pi_m, 'hotpink'):
                mr, ma, mangle = m_sphere = [r + step * pi_r, 
                                             a + step * pi_a/(r * np.sin(angle)), 
                                             angle + step * pi_angle/r] 
                mx, my, mz = sph2cart(ma, np.pi/2 - mangle, mr)
                z, mz = z + height, mz + height  
                ax.plot([x, mx],
                        [y, my],
                        zs = [z, mz], 
                        color = color, 
                        linewidth = 0.3)
                next_starting_points.append(m_sphere)
                
        # Next iterations
        while max(np.array(next_starting_points)[:,0]) < r2:
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
                z, mz = z + height, mz + height  
                ax.plot([x, mx],
                        [y, my],
                        zs = [z, mz], 
                        color = color, 
                        linewidth = 0.3)
                next_starting_points.append(m_sphere)
    plt.show()
michell_3d_spinning(args.point1, args.point2) #user input Radial distribution on bowl 1, Azimutal distribution on bowl 1

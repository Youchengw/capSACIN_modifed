import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
Plot = False
def definePlane(x,y,z):
    coords = np.vstack((x,y,z)).T
    com = (np.mean(x), np.mean(y), np.mean(z))
    coords = coords - com
    
    x=coords[:,0]
    y=coords[:,1]
    z=coords[:,2]
    
    if Plot:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        
        ax.scatter(x, y, z, color='red')
        ax.scatter(0,0,0, color='goldenrod')
    
    ux, uy, uz = u = [x[1]-com[0], y[1]-com[1], z[1]-com[2]]
    vx, vy, vz = v = [x[2]-com[0], y[2]-com[1], z[2]-com[2]]
    # u_cross_v = [uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx]
    u_cross_v = np.cross(u,v)
    
    point = np.array(com)
    normal = np.array(u_cross_v)
    d = -point.dot(normal)
    
    xx, yy = np.meshgrid((-10,10), (-10,10))
    z = (-normal[0]*xx - normal[1]*yy - d) * 1 / normal[2]
    
    if Plot:
        ax.plot_surface(xx,yy,z, alpha=0.6)
    # ax.scatter(normal[0], normal[1], normal[2], color='indigo')
    normalVectorLine = np.vstack((point,normal))
    
    normalVector = np.asarray((point - normal), dtype=float)
    
    if Plot:
        ax.plot(normalVectorLine[:,0], normalVectorLine[:,1], normalVectorLine[:,2], color='indigo')
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
    
        lim=80
        ax.set_xlim(-1*lim,lim)
        ax.set_ylim(-1*lim,lim)
        ax.set_zlim(-1*lim,lim)
        plt.show()
    
    return normal
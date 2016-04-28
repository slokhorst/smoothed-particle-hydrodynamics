# import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import mpl_toolkits.mplot3d
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

class AnimatedScatter(object):
    def __init__(self, pos, domain, updfunc, **updargs):
        self.pos = pos
        self.domain = domain
        self.updfunc = updfunc
        self.updargs = updargs
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ani = ani.FuncAnimation(self.fig, self.update, interval=1, frames=500, init_func=self.setup, blit=False)
        #self.ani.save("sph.mp4", fps=30)

    def setup(self):
        self.scat = self.ax.scatter(self.pos[:,0],self.pos[:,1],self.pos[:,2])
        self.scat._offsets3d = np.transpose(self.pos)
        self.ax.set_xlim(-self.domain, self.domain)
        self.ax.set_ylim(-self.domain, self.domain)
        self.ax.set_zlim(-self.domain, self.domain)
        #sphere
        R = 2
        r0 = np.array([0.0,-self.domain,-self.domain])
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = R*np.outer(np.cos(u), np.sin(v)) + r0[0]
        y = R*np.outer(np.sin(u), np.sin(v)) + r0[1]
        z = R*np.outer(np.ones(np.size(u)), np.cos(v)) + r0[2]
        self.ax.plot_surface(x, y, z, rstride=4, cstride=4, color='r')
        #triangle
        va = np.array([[self.domain,0.0,-self.domain], [-self.domain,0.0,-self.domain], [self.domain,self.domain,0.0]])
        vb = np.array([[-self.domain,self.domain,0], [-self.domain,0.0,-self.domain], [self.domain,self.domain,0.0]])
        vv = np.concatenate([va,vb])
        x = vv[:,0]
        y = vv[:,1]
        z = vv[:,2]
        verts = [list(zip(x, y,z))]
        self.ax.add_collection3d(Poly3DCollection(verts))

        return self.scat,

    def update(self, i):
        self.updfunc(**self.updargs)

        plt.draw()
        return self.scat,

    def show(self):
        plt.show()


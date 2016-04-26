import matplotlib.pyplot as plt
import matplotlib.animation as ani
import mpl_toolkits.mplot3d
import numpy as np

class AnimatedScatter(object):
    def __init__(self, pos, domain, updfunc, **updargs):
        self.pos = pos
        self.domain = domain
        self.updfunc = updfunc
        self.updargs = updargs
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ani = ani.FuncAnimation(self.fig, self.update, interval=1, init_func=self.setup, blit=False)

    def setup(self):
        self.scat = self.ax.scatter(self.pos[:,0],self.pos[:,1],self.pos[:,2])
        self.scat._offsets3d = np.transpose(self.pos)
        self.ax.set_xlim(-self.domain, self.domain)
        self.ax.set_ylim(-self.domain, self.domain)
        self.ax.set_zlim(-self.domain, self.domain)
        return self.scat,

    def update(self, i):
        self.updfunc(**self.updargs)
        plt.draw()
        return self.scat,

    def show(self):
        plt.show()


import numpy as np
import numpy.random
import anim_md

class SPHSim:
	gamma = 7
	h = 1
	N = 100
	rho0 = 1
	C = 0.45
	dt = 0.1
	kernel_const = 315/(64*np.pi*h**9)

	def kernel(self, dr_len):
		return (dr_len < self.h) * self.kernel_const * (self.h**2 - dr_len**2)**3

	def gradkernel(self, dr):
		dr_len = np.sqrt(np.sum(np.square(dr)))
		if dr_len > self.h:
			return 0*dr
		return -6 * dr * self.kernel_const * (self.h**2 - dr_len**2)**2

	def density(self, r):
		dr_len = np.sqrt(np.sum(np.square(r-self.pos),axis=1))
		return np.sum(self.mass*self.kernel(dr_len))

	def pressure(self, rho):
		return self.rho0 * self.C**2 / self.gamma * ((rho/self.rho0)**self.gamma - 1)

	def boundary_force(self, r):
		f = np.array([0.0,0.0,0.0])
		f[2] += -1 #gravity
		f[2] += 10*np.power(r[2],-1) #floor
		return f


	def __init__(self):
		self.pos = np.random.uniform(low=0.0, high=10.0, size=(self.N,3))
		self.vel = np.zeros(shape=(self.N,3))
		self.acc = np.zeros(shape=(self.N,3))
		self.mass = np.ndarray([self.N])
		self.rho = np.ndarray([self.N])
		self.P = np.ndarray([self.N])

		self.mass[:] = 1
		for i in range(0,self.N):
			self.rho[i] = self.density(self.pos[i,:])
		self.P[:] = self.pressure(self.rho)

	def update(self):
		self.vel += 0.5*self.acc*self.dt
		self.pos += self.vel*self.dt

		self.acc[:,:] = 0
		for i in range(0,self.N):
			self.acc[i,:] = 0
			for j in range(0,self.N):
				self.acc[i,:] += -self.mass[j]*((self.P[i]/self.rho[i]**2)+(self.P[j]/self.rho[j]**2))*self.gradkernel(self.pos[i,:]-self.pos[j,:])
			self.acc[i,:] += self.boundary_force(self.pos[i,:])

		self.vel += 0.5*self.acc*self.dt


		for i in range(0,self.N):
			self.rho[i] = self.density(self.pos[i,:])
		self.P[:] = self.pressure(self.rho)


sph = SPHSim()

ascat = anim_md.AnimatedScatter(sph.pos, 100, sph.update)
ascat.show()

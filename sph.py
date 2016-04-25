import numpy as np
import numpy.random
import anim_md

class SPHSim:
	gamma = 7
	h = 1
	N = 100
	rho0 = 1
	C = 0.45
	dt = 0.2
	D = 20
	mu = 1

	kernel_density_const =  315/(64*np.pi*h**9)
	gradkernel_pressure_const = -6 * 315/(64*np.pi*h**9)
	grad2kernel_viscosity_const = 45/(np.pi*h**6)

	def kernel_density(self, dr_len):
		return (dr_len < self.h) * self.kernel_density_const * (self.h**2 - dr_len**2)**3

	def gradkernel_pressure(self, dr_len, dr):
		#return (dr_len < self.h) * -135/(np.pi*self.h**6) * (self.h - dr_len)**2
		return (dr_len < self.h) * dr * self.gradkernel_pressure_const * (self.h**2 - dr_len**2)**2

	def grad2kernel_viscosity(self, dr_len):
		return (dr_len < self.h) * self.grad2kernel_viscosity_const * (self.h - dr_len)

	def density(self, r):
		dr_len = np.sqrt(np.sum(np.square(r-self.pos),axis=1))
		return np.sum(self.mass*self.kernel_density(dr_len))

	def pressure(self, rho):
		return self.rho0 * self.C**2 / self.gamma * ((rho/self.rho0)**self.gamma - 1)

	def body_force(self, r):
		f = np.array([0.0,0.0,0.0])
		#f[2] -= 1 #gravity
		return f


	def __init__(self):
		self.pos = np.random.uniform(low=-10.0, high=10.0, size=(self.N,3))
		self.vel = np.random.uniform(low=-10.0, high=10.0, size=(self.N,3))
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

		self.vel *= 1-2*(np.absolute(self.pos)>self.D)

		self.acc[:,:] = 0
		for i in range(0,self.N):
			self.acc[i,:] = 0
			for j in range(0,self.N):
				dr = self.pos[i,:]-self.pos[j,:]
				dr_len = np.sqrt(np.sum(np.square(dr)))
				if dr_len > self.h:
					continue
				self.acc[i,:] += -self.mass[j]*((self.P[i]/self.rho[i]**2)+(self.P[j]/self.rho[j]**2))*self.gradkernel_pressure(dr_len,dr)
				#self.acc[i,:] += -self.mu*((self.vel[j]/self.rho[j]**2)-(self.vel[i]/self.rho[i]**2))*self.grad2kernel_viscosity(dr_len)
			self.acc[i,:] += self.body_force(self.pos[i,:])

		self.vel += 0.5*self.acc*self.dt


		for i in range(0,self.N):
			self.rho[i] = self.density(self.pos[i,:])
		self.P[:] = self.pressure(self.rho)

		print("total mass:",np.sum(self.rho))
		print("total velocity:",np.sum(np.absolute(self.vel)))


sph = SPHSim()

#ascat = anim_md.AnimatedScatter(sph.pos, 100, sph.update)
#ascat.show()

for i in range(0,100):
	sph.update()

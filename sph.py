import numpy as np
import numpy.random
import anim_md
import os

class SPHSim:
	gamma = 7
	h = 1
	N = 1000
	rho0 = 1
	C = 0.45
	dt = 0.01
	D = 7
	mu = 1
	G = 10
	mass = 1
	interaction_interval = 10

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
		f[2] -= self.G #gravity
		return f

	def update_interacting_particles(self):
	    v_max_sq = 0
	    for i in range(0, self.N):
	    	v_max_sq = max(v_max_sq,np.sum(np.square(self.pos[i,:])))

	    cutoff_sq = self.interaction_interval*self.interaction_interval*self.dt*self.dt*v_max_sq

	    self.interacting = np.zeros((self.N,self.N))
	    for i in range(0, self.N):
	        for j in range(i+1, self.N):
	            dr_sq = np.sum(np.square(self.pos[i,:]-self.pos[j,:]))

	            if dr_sq < cutoff_sq:
	                self.interacting[i,j] = 1
	                self.interacting[j,i] = 1

	def save_pos(self, filename):
		np.savetxt(filename,self.pos,delimiter=',',newline=',\n')

	def __init__(self):
		self.it = 0
		self.pos = np.random.uniform(low=-self.D, high=self.D, size=(self.N,3))
		#self.vel = np.random.uniform(low=-10.0, high=10.0, size=(self.N,3))
		self.vel = np.zeros(shape=(self.N,3))
		self.acc = np.zeros(shape=(self.N,3))
		self.rho = np.ndarray([self.N])
		self.P = np.ndarray([self.N])

		for i in range(0,self.N):
			self.rho[i] = self.density(self.pos[i,:])
		self.P[:] = self.pressure(self.rho)

		self.update_interacting_particles()

	def update(self):
		self.it += 1
		if self.it%self.interaction_interval == 0:
			self.update_interacting_particles()

		self.vel += 0.5*self.acc*self.dt
		self.pos += self.vel*self.dt

		self.vel *= 1-2*(np.absolute(self.pos)>self.D)
		self.pos[:] = np.clip(self.pos,-self.D,self.D)

		self.acc[:,:] = 0
		for i in range(0,self.N):
			for j in np.nonzero(self.interacting[i,:])[0]:
				if j > i:
					dr = self.pos[i,:]-self.pos[j,:]
					dr_len = np.sqrt(np.sum(np.square(dr)))
					if dr_len > self.h:
						continue
					acc_P_i = -self.mass*((self.P[i]/self.rho[i]**2)+(self.P[j]/self.rho[j]**2))*self.gradkernel_pressure(dr_len,dr)
					acc_mu_i = +self.mu*((self.vel[j]/self.rho[j]**2)-(self.vel[i]/self.rho[i]**2))*self.grad2kernel_viscosity(dr_len)
					acc_i = acc_P_i + acc_mu_i
					acc_j = -acc_i
					self.acc[i,:] += acc_i
					self.acc[j,:] += acc_j
			self.acc[i,:] += self.body_force(self.pos[i,:])

		self.vel += 0.5*self.acc*self.dt

		for i in range(0,self.N):
			self.rho[i] = self.density(self.pos[i,:])
		self.P[:] = self.pressure(self.rho)


sph = SPHSim()

# ani = anim_md.AnimatedScatter(sph.pos, sph.D, sph.update)
# ani.show()

os.makedirs('data',exist_ok=True)
for i in range(0,1000):
	print(i)
	sph.update()
	if i%10 == 0:
		sph.save_pos('data/posdata{:03d}.csv'.format(round(i/10)))

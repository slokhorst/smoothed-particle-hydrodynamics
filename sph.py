import numpy as np
import numpy.random
import anim_md
import os
from ray_tri_isec import triangle_intersection, triangle_normal

class Triangle:
	def __init__(self, v1, v2, v3):
		self.v1 = np.array(v1)
		self.v2 = np.array(v2)
		self.v3 = np.array(v3)
		self.n = triangle_normal(self.v1,self.v2,self.v3)

class Sphere:
	def __init__(self, r0, R):
		self.r0 = np.array(r0)
		self.R = R


class SPHSim:
	gamma = 7
	gamma_surf = 1
	h = 1
	nx = 1*8
	nb = 0*2
	N = int(0.25*nx**3) + nb**3
	rho0 = 1
	C = 0.45
	dt = 0.01
	D = 1*3
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

	def normal_vec(self, dr_len, dr):
		return (dr_len < self.h) * dr * self.gradkernel_pressure_const* (self.h**2 - dr_len**2)**2

	def div_normal_vec(self, dr_len):
		return (dr_len < self.h) * self.kernel_density_const*(self.h**2-dr_len**2)*(24*dr_len**2-18*(self.h**2-dr_len**2))

	def density(self, r):
		dr_len = np.sqrt(np.sum(np.square(r-self.pos),axis=1))
		return np.sum(self.mass*self.kernel_density(dr_len))

	def pressure(self, rho):
		return self.rho0 * self.C**2 / self.gamma * ((rho/self.rho0)**self.gamma - 1)

	def body_force(self, r):
		f = np.array([0.0,0.0,0.0])
		f[2] -= self.G #gravity
		return f

	def boundary_interaction(self):
		# make sure the particles stay in the box
		self.vel *= 1-2*(np.absolute(self.pos)>self.D)
		self.pos[:] = np.clip(self.pos,-self.D,self.D)

		# handle collisions with objects
		for i in range(0,self.N):
			for geo in self.geometry:
				if type(geo) == Triangle:
					direction = self.vel[i,:]/np.sqrt(np.sum(np.square(self.vel[i,:])))
					isec, d = triangle_intersection(geo.v1, geo.v2, geo.v3 ,self.pos[i,:],direction)
					if isec and d < 0.1:
						self.vel[i,:] = self.vel[i,:] - 2*(np.sum(self.vel[i,:]*geo.n))*geo.n
				if type(geo) == Sphere:
					dr = self.pos[i,:]-geo.r0
					dr_len = np.sqrt(np.sum(np.square(dr)))
					if dr_len < geo.R:
						n = -dr/dr_len
						self.vel[i,:] = self.vel[i,:] - 2*(np.sum(self.vel[i,:]*n))*n		

	def update_interacting_particles(self):
		v_max_sq = np.max(np.sum(np.square(self.vel),axis=1))
		cutoff_sq = self.interaction_interval**2*self.dt**2*v_max_sq+self.h**2
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
		self.vel = np.zeros(shape=(self.N,3))
		self.pos = np.zeros(shape=(self.N,3))
		i = 0
		for x in np.linspace(-self.D,self.D,num=self.nx):
			for z in np.linspace(0.1,self.D,num=self.nx):
				for y in np.linspace(+self.D,+self.D/2,num=self.nx/4):
					self.pos[i,0] =x
					self.pos[i,1] =y
					self.pos[i,2] =z
					i +=1
		for x in np.linspace(-0.25*self.D,0.25*self.D,num=self.nb):
			for y in np.linspace(-0.25*self.D,0.25*self.D,num=self.nb):
				for z in np.linspace(0.5*self.D,1.0*self.D,num=self.nb):
					self.pos[i,0] =x
					self.pos[i,1] =y
					self.pos[i,2] =z
					i +=1
		self.vel = np.random.uniform(low=-1, high=1, size=(self.N,3))
		self.acc = np.zeros(shape=(self.N,3))
		self.rho = np.ndarray([self.N])
		self.P = np.ndarray([self.N])

		self.geometry = np.array([
			Sphere([0.0,-self.D/2,-self.D],2.0),
			Triangle([self.D,0.0,-self.D], [-self.D,0.0,-self.D], [self.D,self.D,0.0]),
			Triangle([-self.D,self.D,0], [-self.D,0.0,-self.D], [self.D,self.D,0.0]),
		])

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

		self.boundary_interaction()

		self.acc[:,:] = 0
		for i in range(0,self.N):
			for j in np.nonzero(self.interacting[i,:])[0]:
				if j > i:
					dr = self.pos[i,:]-self.pos[j,:]
					dr_len = np.sqrt(np.sum(np.square(dr)))
					if dr_len > self.h:
						continue
					acc_s = 0
					n = self.normal_vec(dr_len,dr)
					len_n = np.sqrt(np.sum(np.square(n)))
					acc_P_i = -self.mass*((self.P[i]/self.rho[i]**2)+(self.P[j]/self.rho[j]**2))*self.gradkernel_pressure(dr_len,dr)
					acc_mu_i = +self.mu*((self.vel[j]/self.rho[j]**2)-(self.vel[i]/self.rho[i]**2))*self.grad2kernel_viscosity(dr_len)
					if len_n > 0:
						acc_s = +self.gamma_surf*(self.div_normal_vec(dr_len)*n/len_n)
					acc_i = acc_P_i + acc_mu_i+acc_s
					acc_j = -acc_i
					self.acc[i,:] += acc_i
					self.acc[j,:] += acc_j
			self.acc[i,:] += self.body_force(self.pos[i,:])

		self.vel += 0.5*self.acc*self.dt

		for i in range(0,self.N):
			self.rho[i] = self.density(self.pos[i,:])
		self.P[:] = self.pressure(self.rho)


sph = SPHSim()

ani = anim_md.AnimatedScatter(sph.pos, sph.D, sph.update)
ani.show()

# os.makedirs('data',exist_ok=True)
# for i in range(0,1000):
# 	print(i)
# 	sph.update()
# 	if i%10 == 0:
# 		sph.save_pos('data/posdata{:03d}.csv'.format(round(i/10)))

import numpy as np
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
	sigma = 1
	h = 1
	rho0 = 1
	C = 0.45
	dt = 0.01
	mu = 1
	G = 10
	mass = 1
	damping = 0.0
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
		self.vel *= 1-(2-self.damping)*(np.absolute(self.pos)>self.D)
		self.pos[:] = np.clip(self.pos,-self.D,self.D)

		# handle collisions with objects
		for i in range(0,self.N):
			for geo in self.geometry:
				if type(geo) == Triangle:
					direction = self.vel[i,:]/np.sqrt(np.sum(np.square(self.vel[i,:])))
					isec, d = triangle_intersection(geo.v1, geo.v2, geo.v3 ,self.pos[i,:],direction)
					ds_sq = np.sum(np.square(self.vel[i,:]*self.dt + 0.5*self.acc[i,:]*self.dt**2))
					if isec and (d)**2 < 2*ds_sq:
						self.vel[i,:] = self.vel[i,:] - (2-self.damping)*(np.sum(self.vel[i,:]*geo.n))*geo.n
				if type(geo) == Sphere:
					dr = self.pos[i,:]-geo.r0
					dr_len = np.sqrt(np.sum(np.square(dr)))
					if dr_len < geo.R:
						n = -dr/dr_len
						self.vel[i,:] = self.vel[i,:] - (2-self.damping)*(np.sum(self.vel[i,:]*n))*n		

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
		np.savetxt(filename,self.pos/self.D,delimiter=',',newline=',\n')

	def __init__(self, N, D):
		self.N = N
		self.D = D
		self.it = 0
		self.vel = np.zeros(shape=(self.N,3))
		self.pos = np.zeros(shape=(self.N,3))
		self.acc = np.zeros(shape=(self.N,3))
		self.rho = np.ndarray([self.N])
		self.P = np.ndarray([self.N])
		self.geometry = []

	def update(self):
		if self.it%self.interaction_interval == 0:
			self.update_interacting_particles()

		for i in range(0,self.N):
			self.rho[i] = self.density(self.pos[i,:])
		self.P[:] = self.pressure(self.rho)

		self.vel += 0.5*self.acc*self.dt
		self.pos += self.vel*self.dt

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
						acc_s = +self.sigma*(self.div_normal_vec(dr_len)*n/len_n)
					acc_i = acc_P_i + acc_mu_i+acc_s
					acc_j = -acc_i
					self.acc[i,:] += acc_i
					self.acc[j,:] += acc_j
			self.acc[i,:] += self.body_force(self.pos[i,:])

		self.vel += 0.5*self.acc*self.dt

		self.boundary_interaction()

		self.it += 1

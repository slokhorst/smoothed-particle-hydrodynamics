import numpy as np
import numpy.random
import anim_md
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class SPHSim:
	gamma = 7
	gamma_surf = 4
	h = 2
	h_s = 1
	N = 100
	rho0 = 1
	C = 0.45
	dt = 0.1
	D = 4
	mu = 40
	G = 1
	mass = 1
	interaction_interval = 10
	sigma = 10

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

	def adh_kernel(self, dr_len):
		
		y=0

		if dr_len < 0.5*self.h_s:
			y = 32/(np.pi*self.h_s**9)*(2*(self.h_s-dr_len)**3*dr_len**3-self.h_s**6/64)
		else:
			y=  32/(np.pi*self.h_s**9)*(self.h_s-dr_len)**3*dr_len**3

		return y

	def normal_vec(self, dr_len, dr):
		return (dr_len < self.h) * dr * self.gradkernel_pressure_const* (self.h**2 - dr_len**2)**2

	def div_normal_vec(self, dr_len):
		return (dr_len < self.h) * self.kernel_density_const*(self.h**2-dr_len**2)*(24*dr_len**2-18*(self.h**2-dr_len**2))

	def init_pos_box(self):
		W = 5
		L = 5
		Dp = 4

		horiz_shift = -0.4*self.D*np.array([1,1,0])
		vert_shift = 0.3*self.D*np.array([0,0,1])

		init_pos = np.zeros((L*W*Dp,3))

		offs = 0
		for i in range(0,Dp):
			for j in range(0,W):
				for k in range(0,L):
					init_pos[offs+i+j+k,:] = horiz_shift + vert_shift + k*1/6*self.D*np.array([1,0,0])\
					+ j*1/6*self.D*np.array([0,1,0])- i*1/5*self.D*np.array([0,0,1])
				offs += L-1
			offs += W-1

		return init_pos

	def density(self, r):
		dr_len = np.sqrt(np.sum(np.square(r-self.pos),axis=1))
		return np.sum(self.mass*self.kernel_density(dr_len))

	def pressure(self, rho):
		return self.rho0 * self.C**2 / self.gamma * ((rho/self.rho0)**self.gamma - 1)

	def body_force(self, r):
		f = np.array([0.0,0.0,0.0])
		#f[2] -= self.G #gravity
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
		self.pos = np.zeros((self.N,3))
		# self.pos[:,0:2] = np.random.uniform(low=-self.D, high=self.D, size=(self.N,2))
		# self.pos[:,2] = np.random.uniform(low=-self.D, high=-0.7*self.D, size=(self.N))
		self.pos = self.init_pos_box()
		#self.vel = np.random.uniform(low=-1.0, high=1.0, size=(self.N,3))
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
		self.pos = np.clip(self.pos,-self.D,self.D)
		self.acc[:,:] = 0

		for i in range(0,self.N):
			for j in range(0,self.N):#np.nonzero(self.interacting[i,:])[0]:
				# if j > i:
				dr = self.pos[i,:]-self.pos[j,:]
				dr_len = np.sqrt(np.sum(np.square(dr)))
				#if dr_len > self.h and dr_len > 2*self.h_s:
				# 	continue
				acc_s = 0
				n = self.normal_vec(dr_len,dr)
				len_n = np.sqrt(np.sum(np.square(n)))
				acc_P_i = -self.mass*((self.P[i]/self.rho[i]**2)+(self.P[j]/self.rho[j]**2))*self.gradkernel_pressure(dr_len,dr)
				acc_mu_i = +self.mu*((self.vel[j]/self.rho[j]**2)-(self.vel[i]/self.rho[i]**2))*self.grad2kernel_viscosity(dr_len)
				if len_n > 0.2:
					acc_s = +self.gamma_surf*(self.div_normal_vec(dr_len)*n/len_n)

				#acc_coh = -self.gamma_surf*self.adh_kernel(dr_len)*(dr/dr_len)

				# print('P_i={0} mu={1} coh={2}'.format(acc_P_i,acc_mu_i,acc_coh))

				acc_i = (acc_P_i + acc_mu_i + acc_s)
				# acc_j = -acc_i
				self.acc[i,:] += acc_i
				# self.acc[j,:] += acc_j
			
			self.acc[i,:] += self.body_force(self.pos[i,:])
		print('total velocity=',np.sum(np.absolute(self.vel)))
		self.vel += 0.5*self.acc*self.dt

		for i in range(0,self.N):
			self.rho[i] = self.density(self.pos[i,:])
		self.P[:] = self.pressure(self.rho)


sph = SPHSim()

ascat = anim_md.AnimatedScatter(sph.pos, 4, sph.update)
ascat.show()

# for i in range(0,10000):
# 	print(i)
# 	fig = plt.figure()
# 	ax = fig.add_subplot(111, projection='3d')
# 	ax.scatter(sph.pos[:,0], sph.pos[:,1], sph.pos[:,2])
# 	plt.show()
# 	sph.update()

# sph.save_pos('posdata')
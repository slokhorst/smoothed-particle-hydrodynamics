import anim_md
import os
import numpy as np
import numpy.random
from sph import SPHSim, Triangle, Sphere

nx = 0*8
nb = 2*2
N = int(0.25*nx**3) + nb**3
D = 70

sph = SPHSim(N,D)

sph.geometry.append(Triangle([D,D,-10], [-D,10,-D], [D,10,-D]))
sph.geometry.append(Triangle([-D,D,-10], [-D,10,-D], [D,D,-10]))
sph.geometry.append(Triangle([-5,-D,-D], [-5,D,D], [-5,D,-D]))
sph.geometry.append(Triangle([-5,-D,-D], [-5,-D,D], [-5,D,D]))
sph.geometry.append(Triangle([5,-D,-D], [5,D,D], [5,D,-D]))
sph.geometry.append(Triangle([5,-D,-D], [5,-D,D], [5,D,D]))

sph.geometry.append(Triangle([5,-20,-D], [5,-30,-D+10], [-5,-20,-D]))
sph.geometry.append(Triangle([-5,-20,-D], [-5,-30,-D+10], [5,-20,-D]))

i = 0
for x in np.linspace(-D,D,num=nx):
	for z in np.linspace(0.1,D,num=nx):
		for y in np.linspace(+D,+D/2,num=nx/4):
			sph.pos[i,0] =x
			sph.pos[i,1] =y
			sph.pos[i,2] =z
			i +=1
for x in np.linspace(-4,4,num=nb):
			for y in np.linspace(D-4,D ,num=nb):
				for z in np.linspace(-10,2*nb-10,num=nb):
					sph.pos[i,0] =x
					sph.pos[i,1] =y
					sph.pos[i,2] =z
					i +=1
sph.vel = np.random.uniform(low=-1, high=1, size=(N,3))

ani = anim_md.AnimatedScatter(sph.pos, D, sph.update)

for geo in sph.geometry:
	if type(geo) == Sphere:
		ani.add_sphere(geo)
	if type(geo) == Triangle:
		ani.add_triangle(geo)

ani.show()

# os.makedirs('data',exist_ok=True)
# for i in range(0,1000):
# 	print(i)
# 	sph.update()
# 	if i%10 == 0:
# 		sph.save_pos('data/posdata{:03d}.csv'.format(round(i/10)))
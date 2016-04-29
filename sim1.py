import anim_md
import os
import numpy as np
import numpy.random
from sph import SPHSim, Triangle, Sphere

nx = 12
N = int(nx**3/2)
D = 10*1

sph = SPHSim(N,D)
sph.mu = 10
sph.damping = 0.5

sph.geometry.append(Sphere([0.0,-D/2,-D],0.3*D))
sph.geometry.append(Triangle([D,D,0.0], [-D,0.0,-D], [D,0.0,-D]))
sph.geometry.append(Triangle([-D,D,0], [-D,0.0,-D], [D,D,0.0]))

i = 0
for x in np.linspace(-D,D,num=2*nx):
	for z in np.linspace(3/4*D,1/4*D,num=nx/2):
		for y in np.linspace(+D,+D/2,num=nx/2):
			sph.pos[i,0] =x
			sph.pos[i,1] =y
			sph.pos[i,2] =z
			i +=1
sph.vel = np.random.uniform(low=-0.1, high=0.1, size=(N,3))

ani = anim_md.AnimatedScatter(sph.pos, D, sph.update)
for geo in sph.geometry:
	if type(geo) == Sphere:
		ani.add_sphere(geo)
	if type(geo) == Triangle:
		ani.add_triangle(geo)
#ani.save("sim1.mp4")
ani.show()

# os.makedirs('data',exist_ok=True)
# for i in range(0,1000):
# 	print(i)
# 	sph.update()
# 	if i%10 == 0:
# 		sph.save_pos('data/posdata{:03d}.csv'.format(round(i/10)))

import numpy as np
EPSILON=0.000001

def cross(v1, v2):
	return np.array([v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0]])
def dot(v1, v2):
	return np.sum(v1*v2)

def triangle_normal(v1, v2, v3):
	e1 = v2-v1
	e2 = v3-v1
	normal = cross(e1,e2)
	return normal/np.sqrt(np.sum(np.square(normal)))

def triangle_intersection(v1, v2, v3, pos, vel):
	# Find vectors for two edges sharing V1
	e1 = v2-v1
	e2 = v3-v1
	# Begin calculating determinant - also used to calculate u parameter
	P = cross(vel, e2)
	# if determinant is near zero, ray lies in plane of triangle
	det = dot(e1, P)
	# NOT CULLING
	if(det > -EPSILON and det < EPSILON):
		return False,0
	inv_det = 1.0 / det

	# calculate distance from V1 to ray origin
	T = pos - v1

	# Calculate u parameter and test bound
	u = dot(T, P) * inv_det
	# The intersection lies outside of the triangle
	if(u < 0.0 or u > 1.0):
		return False,0

	# Prepare to test v parameter
	Q = cross(T, e1)

	# Calculate V parameter and test bound
	v = dot(vel, Q) * inv_det
	# The intersection lies outside of the triangle
	if(v < 0.0 or u + v  > 1.0):
		return False,0

	t = dot(e2, Q) * inv_det

	if(t > EPSILON): # ray intersection
		return True,t

	#  No hit, no win
	return False,0
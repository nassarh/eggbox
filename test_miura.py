# coding=UTF8

import numpy as np
import miura as mr

# Pour matplotlib
# from mpl_toolkits.mplot3d import axes3d
# import matplotlib.pyplot as plt
# from matplotlib import colors
# import matplotlib.cm as cm

# Pour mayavi
import mayavi.mlab as mayalab

###############################################
###############################################
e1 = np.array([1.,0,0])
e2 = np.array([0,1.,0])
e3 = np.array([0,0,1.])

cas_test = 1

if cas_test == 0:
	# propager une droite
	# l'instabilité se manifeste à partir de n=170
	n=10
	r=1
	theta=0
	phi=-np.pi/2*0.995

	a,u,v = mr.zigzag(phi,n=400)
else: 
	# propager un cercle
	r=1.
	rayon = 100.
	theta = np.pi/3
	a,u,v = mr.zigcercle(theta,rayon)
	#a,u,v = mr.twist(theta,np.pi/30)

X, Y, Z = mr.xyz_emission(a,u,v)

#X = X - X.mean()
#Y = Y - Y.mean()
#Z = Z - Z.mean()


#X = X[::4,::4]
#Y = Y[::4,::4]
#Z = Z[::4,::4]

x = X.ravel()
y = Y.ravel()
z = Z.ravel()

########################################################################################
# Utiliser
# bloc 1 pour visualiser avec mayavi
# bloc 2 pour visualiser avec matplotlib
# bloc 3 pour exporter un fichier .vtk (pour paraview)
########################################################################################
##########################
# BLOC 1##################
##########################
mesh = mayalab.triangular_mesh(x,y,z,mr.triangles(X.shape[0]),color=(0,116./255,64./255))

mesh.actor.property.interpolation = 'phong'
mesh.actor.property.specular = 1
mesh.actor.property.specular_power = 128
mesh.actor.property.edge_visibility = 1
mesh.actor.property.ambient = 1

mayalab.show()

##########################
# BLOC 2##################
##########################
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_trisurf(x,y,np.array(mr.triangles(X.shape[0])),z,shade=True,color=(0.2,0.4,1))

# m = min(np.min(X),np.min(Y),np.min(Z))
# M = max(np.max(X),np.max(Y),np.max(Z))
# ax.set_zlim(m,M)
# plt.xlim(m,M)
# plt.ylim(m,M)
# ax.axis('off')
# plt.axis('equal')

# plt.show()
##########################
# BLOC 3##################
##########################
# P = np.vstack((x,y,z)).T

# mr.writeVTK("test.vtk",P,mr.triangles(X.shape[0]))

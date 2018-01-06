# coding=UTF8

import numpy as np
import numpy.linalg as npal

#########################################
# Construire l'intersection de 3 sphères
#########################################
def wedge(u,v):
	# retourne le produit vectoriel de u et v
	return np.cross(u,v)

def w(u,v,sgn=1):
        # construit l'intersection de trois sphères unités centrées sur 0, u et v
        # sgn selectionne l'une des deux solutions possibles
	#u = u/npal.norm(u)
	#v = v/npal.norm(v)
	C = np.dot(u,v)
	c = 1./np.sqrt(2*(1+C))
	# solution existe sauf
        if C < -0.5 :
		raise ValueError('Pas de solution')
	else:
                s = np.sqrt((1+2*C)/2/(1+C))
	
        w = c*(u+v)/npal.norm(u+v) + sgn*s*wedge(u,v)/npal.norm(wedge(u,v))
	#w = c*(u+v)/np.sqrt(2+2*C) + sgn*s*wedge(u,v)/npal.norm(wedge(u,v))
	
	#w = w/npal.norm(w)
	
	return w
	
#########################################
# Problème "de Cauchy"
#########################################
def emission_etape(a,u,v,sgn=1):
        # itérer la construction 3-sphères pour un zigzag a
	N = np.shape(a)[0]-2
	aa = np.copy(a[1:-1,:])
	ww = np.copy(u)
        if N > 1:
            uu = np.copy(u[1:,:])
            vv = np.copy(v[1:,:])
        else:
            uu = np.zeros((2,3))
            vv = np.zeros((2,3))
        if N < 1:
            return aa, uu, vv

        p = a[1,:]
		
	wwtemp = w(u[0,:],v[0,:],sgn=sgn)
	ww[0,:] = wwtemp
	
	aa[2*0,:] = p + wwtemp
	
        vv[0,:] = -u[0,:] + wwtemp

        for n in range(1,(N-1)//2):
		p = a[2*n+1,:]
		
		wwtemp = w(u[n,:],v[n,:],sgn=sgn)
		ww[n,:] = wwtemp
		
		aa[2*n,:] = p + wwtemp
		
                uu[n-1,:] = -v[n,:] + wwtemp
                vv[n,:] = -u[n,:] + wwtemp

        n = (N-1)//2
        p = a[2*n+1,:]		
        wwtemp = w(u[n,:],v[n,:],sgn=sgn)
        ww[n,:] = wwtemp

        aa[2*n,:] = p + wwtemp
        uu[n-1,:] = -v[n,:] + wwtemp

	return aa, uu, vv

def emission(a,u,v):
        # intérer la construction des zigzag; a est le zigzag initial
	N = (np.shape(a)[0]+1)//2
	phi = np.zeros((N,3*N))
	aa = np.copy(a)
        uu = np.copy(u)
        vv = np.copy(v)
	
	nosing = 1
	
	for n in range(0,N):
		phi[0:N-n,3*n:3*n+3] = aa[::2,:]
                #print("n="+str(n)+" z="+str(phi[0,3*n+2]))
		
		if nosing :
			try :
				aa, uu, vv = emission_etape(aa,uu,vv,sgn=(((n+1)%2)*(-1)**(n//2) + (n%2)*(-1)**((n+1)//2)))
			# si une exception "Pas de solution", on arrête la construction
			except ValueError as e:
				print(e)
				print("apres " + str(n) + " iterations")
				aa = aa[1:-1,:]
				nosing = 0
		else :
			aa = aa[1:-1,:]
	
	return phi

def xyz_emission(a,u,v,r=1):
	# Cette fonction retourne les coordonnees le parametrage (X, Y, Z) de la surface
	# a est une ligne en zigzag
	# r est la distance entre deux points successifs
	
	# on propage a dans une direction
	phi = emission(a,u,v)

	# on propage a dans l'autre direction ("back propagation")
	# il faut alors lire a a l'envers et oublier ses extremites
	aa = a[1:-1,:][::-1,:]
        uu = -u[0:-1,:][::-1,:]
        vv = -v[1:,:][::-1,:]
	phib = emission(aa,uu,vv)
	
	# phib est lu a l'envers, pour le coller a phi, il faut le redresser d'abord
	# on le lit a l'envers dans les deux dimensions
	phib = phib[::-1,::-1]
	
	# Ceci a aussi l'effet d'interchanger les X et les Z, il faut corriger ca
	for n in range(0,np.shape(phib)[1]//3):
			temp = phib[:,3*n].copy()
			phib[:,3*n] = phib[:,3*n+2].copy()
			phib[:,3*n+2] = temp
	
	# On concatene les deux moities
	phi[1:,3:] = phi[1:,3:] + phib
	
	# ordonner les deux moities pour commencer au commencement
	for r in range(0,np.shape(phi)[0]):
		phi[r,:] = np.roll(phi[r,:],3*r)
	
	return phi[:,0::3], phi[:,1::3], phi[:,2::3]

#########################################
# Des outils de dessin
#########################################
def zigzag(phi,r=1.,n=100,x0=0.,y0=0.):
	a = np.zeros((n+1,3))
	
	a[:,0] = x0+np.linspace(0,n,n+1)*np.cos(phi)
	a[::2,1] = -np.sin(phi)
        a[:,1] = a[:,1]+y0
        #print(a)

        u = np.zeros((n//2,3))
        v = np.zeros((n//2,3))
        #v = a[0:-2:2,:] - a[1:-1:2,:]
        v[:,0] = -np.cos(phi)
        v[:,1] = -np.sin(phi)
        #u = a[2::2,:] - a[1:-1:2,:]
        u[:,0] = np.cos(phi)
        u[:,1] = -np.sin(phi)
	
	return a,u,v

def zigcercle(theta,rayon,r=1.):
	s = np.sin(theta/2)
	c = np.cos(theta/2)
	u0 = np.array([s,c,0])
	v0 = np.array([-s,c,0])
	
	C = 1-2*(s/rayon)**2
	phi0 = np.arccos(C)
	#phi = phi0
	S = np.sin(phi0)
	rot = np.array([[C,0.,-S],[0.,1.,0.],[S,0.,C]])
	
	uv = np.vstack((v0,u0))
	a = np.vstack((np.array([0,rayon,0])+v0,np.array([0,rayon,0]),np.array([0,rayon,0])+u0))
	n = 0
	N = np.floor(np.pi*rayon/2/np.sin(theta))
	#N = 10
	
	while n <= N:
	    #utemp = np.dot(u0,rot)
	    #vtemp = np.dot(v0,rot)
	    utemp = np.dot(uv[-1,:],rot)
	    vtemp = np.dot(uv[-2,:],rot)
	    uv = np.vstack((uv,vtemp,utemp))
	    atemp = a[-1,:]-vtemp
	    a = np.vstack((a,atemp,atemp+utemp))
	    n += 1
	    #phi += phi0
	    #C = np.cos(phi)
	    #S = np.sin(phi)
            #rot = np.array([[C,0.,-S],[0.,1.,0.],[S,0.,C]])

	asym = np.copy(a[3:,:])
	uvsym = np.copy(uv[2:,:])
	asym[:,0] = -asym[:,0]
	uvsym[:,0] = -uvsym[:,0]
	
	a = np.vstack((asym[::-1,:],a))
	uv = np.vstack((uvsym[::-1,:],uv))
	
	return a, uv[1::2,:], uv[::2,:]

def twist(theta,phi0,r=1.):
	s = np.sin(theta/2)
	c = np.cos(theta/2)
	u0 = np.array([s,c,0])
	v0 = np.array([-s,c,0])
	
	phi = phi0
        C = np.cos(phi0)
	S = np.sin(phi0)
	rot = np.array([[1,0,0],[0,C,-S],[0,S,C]])
	
	uv = np.vstack((v0,u0))
	a = np.vstack((v0,np.array([0,0,0]),u0))
	n = 0
	N = 30
	
	while n <= N:
	    utemp = np.dot(u0,rot)
	    vtemp = np.dot(v0,rot)
	    uv = np.vstack((uv,vtemp,utemp))
	    atemp = a[-1,:]-vtemp
	    a = np.vstack((a,atemp,atemp+utemp))
	    n += 1
	    phi += phi0
	    C = np.cos(phi)
	    S = np.sin(phi)
       	    rot = np.array([[1,0,0],[0,C,-S],[0,S,C]])
	
	return a, uv[1::2,:], uv[::2,:]
#########################################
# Des outils de post-traitement
#########################################
def triangles(n):
	return [[i+n*j,i+1+n*j,i+n*(j+1)] for i in range(n-1) for j in range(n-1)] + [[i+n*(j+1),i+1+n*j,i+1+n*(j+1)] for i in range(n-1) for j in range(n-1)]

def lesEs(X,Y,Z):
	# retourne la distance entre les deux extremites des V de tous les zigzags (= 4s^2)
	Ux,Uy,Uz = X[1:,1:] - X[:-1,:-1], Y[1:,1:] - Y[:-1,:-1], Z[1:,1:] - Z[:-1,:-1]
	U = Ux**2+Uy**2+Uz**2
	return U

# test de la fonction appuyer_sur_bord
#if __name__ == "__main__":
    # 4 tests avec une ouverture décroissante
	#emission(zigzag(np.pi/6*1.1,n=2000))
	#emission(zigzag(np.pi/6*1.2,n=2000))
	#emission(zigzag(np.pi/6*1.3,n=2000))
	#emission(zigzag(np.pi/6*1.4,n=2000))
from math import *
import numpy as np

# Geometric parameters
l = 0.2828427
d = 0.4
# Dynamic parameters
mp = 3.0
mf = 1.0


def dgm(q11, q21, assembly_mode):
    x = 0
    y = 0
    q1=q11
    q2=q21
    A2H=np.array([-d/2,(q1-q2)/2])
    
    m = np.array([[0 ,-1], [1, 0]])
    a = np.linalg.norm(A2H)
    print("aaa",a)
    h = sqrt(pow(l,2)-pow(a,2))
    HC= assembly_mode*h/a*m @ A2H
    OA02 = np.array([d/2, 0])
    A02A2= np.array([0, q2])
    OC= OA02+ A02A2 + A2H + HC
    x = OC[0] 
    y = OC[1]
    return x, y


def igm(x, y, gamma1, gamma2):
	q1 = 0
	q2 = 0
	
	q1=y+gamma1*sqrt((l*l)-(x+(d/2))**2)
	q2=y+gamma2*sqrt((l*l)-pow(x-(d/2),2))
	q11=q1
	q21=q2
	return q11, q21


def dgm_passive(q1, q2, assembly_mode):
	q12 = 0.0
	q22 = 0.0
	x,y=dgm(q1, q2, assembly_mode)
	
	phi1=atan2((y-q1)/l,(d/2+x)/l)
	phi2=pi-atan2((y-q2)/l,(d/2-x)/l)
	
	q12=phi1
	q22=phi2
	
	return q12, q22


# You can create intermediate functions to avoid redundant code
def compute_A_B(q11, q12, q21, q22, is_actuated):
	A = 0
	B = 0
	
	q1=q11
	q2=q21
	phi1=q12
	phi2=q22
	
	y0=np.array([0,1])
	U1=np.array([cos(phi1),sin(phi1)])
	U2=np.array([cos(phi2),sin(phi2)])
	
	R=np.array([[0,-1],[1,0]])
	
	v11=R@y0
	v1=R@U1
	v2=R@U2
	
	if is_actuated:
		A=np.array([U1,U2])
		B=np.array([[np.dot(U1,y0),0],[0,np.dot(U2,y0)]])
	else:
		A=np.array([v1,v2])
		B=np.array([[np.dot(v1,y0),0],[0,np.dot(v2,y0)]])
	
	return A, B


def dkm(q11, q12, q21, q22, q11D, q21D):
	xD = 0
	yD = 0
	
	q1=q11
	q2=q21
	phi1=q12
	phi2=q22
	
	A,B=compute_A_B(q1,phi1,q2,phi2,True)
	Ainv=np.linalg.inv(A)
	qD = np.array([q11D, q21D])
	xiD=Ainv @ B @ qD 
	xD=xiD[0]
	yD=xiD[1]
	
	return xD, yD


def ikm(q11, q12, q21, q22, xD, yD):
	q11D = 0
	q21D = 0
	
	A,B=compute_A_B(q11,q12,q21,q22,True)
	Binv=np.linalg.inv(B)
	xiD=np.array([xD,yD])
	qD=Binv @ A @ xiD
	q11D=qD[0]
	q21D=qD[1]
	
	return q11D, q21D


def dkm_passive(q11, q12, q21, q22, q11D, q21D, xD, yD):
	q12D = 0
	q22D = 0
	A,B=compute_A_B(q11,q12,q21,q22,False)
	xiD=np.array([xD,yD])
	qD = np.array([q11D, q21D])
	qdD=(A @ xiD - B @ qD)/l
	q12D= qdD[0]
	q22D= qdD[1]
	return q12D, q22D


def dkm2(q11, q12, q21, q22, q11D, q12D, q21D, q22D, q11DD, q21DD):
	xDD = 0
	yDD = 0
	
	qDD = np.array([q11DD, q21DD])
	A,B=compute_A_B(q11,q12,q21,q22,True)
	Ainv=np.linalg.inv(A)
	d=-l*np.array([q12D*q12D,q22D*q22D])
	xiDD=Ainv @ (B @ qDD + d)
	xDD=xiDD[0]
	yDD=xiDD[1]
	
	return xDD, yDD


def ikm2(q11, q12, q21, q22, q11D, q12D, q21D, q22D, xDD, yDD):
	q11DD = 0
	q21DD = 0
	
	A,B=compute_A_B(q11,q12,q21,q22,True)
	d=-l*np.array([q12D*q12D,q22D*q22D])
	Binv=np.linalg.inv(B)
	xiDD=np.array([xDD,yDD])
	qDD=Binv @ ( A @ xiDD - d)
	q11DD=qDD[0]
	q21DD=qDD[1]
	
	return q11DD, q21DD


def dkm2_passive(q11, q12, q21, q22, q11D, q12D, q21D, q22D, q11DD, q21DD, xDD, yDD):
	q12DD = 0
	q22DD = 0
	
	A,B=compute_A_B(q11,q12,q21,q22,False)
	xiDD=np.array([xDD,yDD])
	qDD = np.array([q11DD, q21DD])
	qdDD=(A @ xiDD - B @ qDD)/l
	q12DD=qdDD[0]
	q22DD=qdDD[1]
	return q12DD, q22DD




def dynamic_model(q11, q12, q21, q22, q11D, q12D, q21D, q22D):
    # Step 1: Compute the mass matrix M1 for the frame
    M1 = np.array([[mf, 0.0], [0.0, mf]])
    
    # Step 2: Compute matrices A and B using the joint angles q11, q12, q21, q22
    A, B = compute_A_B(q11, q12, q21, q22,True)
    
    # Step 3: Invert matrix A to compute Ainv
    Ainv = np.linalg.inv(A)
    
    # Step 4: Compute the Jacobian matrix J from Ainv and B
    J = Ainv@B
    
    # Step 5: Compute the end-effector velocities [xD, yD] using dkm function
    xD, yD = dkm(q11, q12, q21, q22, q11D, q21D)
    
    # Step 6: Form the velocity vector of the end-effector
    xiD = np.array([xD, yD])
    
    # Step 7: Form the velocity vector for the actuated joints
    qaD = np.array([q11D, q21D])
    
    # Step 8: Compute the matrix AD for velocity contributions of the end-effector
    AD = np.array([[2 * xD, 2 * (yD - q11D)],[2 * xD, 2 * (yD - q21D)]])
    
    # Step 9: Compute the matrix BD for velocity contributions of the joints
    BD = np.array([[-2 * (yD - q11D), 0],[0, -2 * (yD - q21D)]])
    
    # Step 10: Compute the vector K using AD, BD, posD, and qaD
    K = AD@xiD + BD@qaD
    
    # Step 11: Compute the vector b as the product of -Ainv and K
    b = -Ainv@K
    
    # Step 12: Compute the full mass matrix M using M1 and the Jacobian matrix J
    M = M1 + mp * np.transpose(J)@ J
    
    # Step 13: Compute the Coriolis vector c using J, Ainv, and b
    c = mp * np.transpose(J)@ b
    
    # Return the mass matrix M and the Coriolis vector c
    return M, c

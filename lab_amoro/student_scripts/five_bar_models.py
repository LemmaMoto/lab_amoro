from math import *
import numpy as np

# Geometric parameters
l = 0.09
d = 0.118
# Dynamic parameters
ZZ1R = 1.0 * 0.045 * 0.045
ZZ2R = 1.0 * 0.045 * 0.045
mR = 0.5


def dgm(q11, q21, assembly_mode):
	x = 0
	y = 0
	A22H=1/2*np.array([-l*cos(q21)-d+l*cos(q11),-l*sin(q21)+l*sin(q11)])
	
	a=np.linalg.norm(A22H)
	
	h=sqrt(l*l-a*a)
	
	HA13=assembly_mode * h/a * np.array([[0,-1],[1,0]]) @ A22H
	
	OA13=np.array([d/2,0])+np.array([l*cos(q21),l*sin(q21)])+A22H+HA13
	
	x=OA13[0]
	y=OA13[1]
	
	return x, y


def igm(x, y, gamma1, gamma2):
    A11A13=np.array([d/2,0])+np.array([x,y])
    A11M1=1/2*A11A13
    c=np.linalg.norm(A11M1) #modulus of A11M1
    b=sqrt(l*l-c*c)
    
    M1A12=gamma1*(b/c)*np.array([[0,-1],[1,0]]) @ A11M1
    
    A11A12=A11M1+M1A12
    
    q11=atan2(A11A12[1]-0,A11A12[0])
    
    #secondary branch
    A21A13=np.array([-d/2,0])+np.array([x,y])
    A21M2=1/2*A21A13
    c=np.linalg.norm(A21M2)
    b=sqrt(l*l-c*c)
    
    M2A22=gamma2*(b/c)*np.array([[0,-1],[1,0]]) @ A21M2
    
    A21A22=A21M2+M2A22
    
    q21=atan2(A21A22[1],A21A22[0])
    return q11, q21


def dgm_passive(q11, q21, assembly_mode):
    q12 = 0.0
    q22 = 0.0
    x,y=dgm(q11, q21, assembly_mode)
    q12= atan2(y/l-sin(q11),x/l + d/(2*l) - cos(q11)) - q11
    q22= atan2(y/l-sin(q21),x/l - d/(2*l) - cos(q21)) - q21
    return q12, q22
	

# You can create intermediate functions to avoid redundant code
def compute_A_B(q11, q12, q21, q22, is_actuated ):
	A = 0
	B = 0
	U11=np.array([cos(q11),sin(q11)])
	U12=np.array([cos(q12+q11),sin(q12+q11)])
	U21=np.array([cos(q21),sin(q21)])
	U22=np.array([cos(q22+q21),sin(q22+q21)])
	#v vectors
	v11=np.array([[0,-1],[1,0]])@U11
	v12=np.array([[0,-1],[1,0]])@U12
	v21=np.array([[0,-1],[1,0]])@U21
	v22=np.array([[0,-1],[1,0]])@U22
	if is_actuated:
		A=np.array([U12,U22])
		B=np.array([[l*np.dot(U12,v11),0],[0,l*np.dot(U22,v21)]])
	else:
		A=np.array([v12,v22])
		B=np.array([[l*np.dot(v12,v11)+l,0],[0,l*np.dot(v22,v21)+l]])
	return A, B
	
def compute_d(q11, q12, q21, q22,q11D, q12D, q21D, q22D, isActuated):
	U11=np.array([cos(q11),sin(q11)])
	U12=np.array([cos(q12+q11),sin(q12+q11)])
	U21=np.array([cos(q21),sin(q21)])
	U22=np.array([cos(q22+q21),sin(q22+q21)])
	#v vectors
	v11=np.array([[0,-1],[1,0]])@U11
	v12=np.array([[0,-1],[1,0]])@U12
	v21=np.array([[0,-1],[1,0]])@U21
	v22=np.array([[0,-1],[1,0]])@U22
	if isActuated:
		d=np.array([-l*pow(q11D,2)*np.dot(U12,U11)-l*pow(q11D+q12D,2),-l*pow(q21D,2)*np.dot(U22,U21)-l*pow(q21D+q22D,2)])
	else:
		d=np.array([-l*pow(q11D,2)*np.dot(v12,U11),-l*pow(q21D,2)*np.dot(v22,U21)])
	return d
	


def dkm(q11, q12, q21, q22, q11D, q21D):
    xD = 0
    yD = 0
    A,B=compute_A_B(q11,q12,q21,q22,True)
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
    d=compute_d(q11, q12, q21, q22,q11D, q12D, q21D, q22D, True)
    xiDD=Ainv @ (B @ qDD + d)
    xDD=xiDD[0]
    yDD=xiDD[1]
    
    return xDD, yDD


def ikm2(q11, q12, q21, q22, q11D, q12D, q21D, q22D, xDD, yDD):
    q11DD = 0
    q21DD = 0
    A,B=compute_A_B(q11,q12,q21,q22,True)
    d=compute_d(q11, q12, q21, q22,q11D, q12D, q21D, q22D, True)
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
    d=compute_d(q11, q12, q21, q22,q11D, q12D, q21D, q22D, False)
    xiDD=np.array([xDD,yDD])
    qDD = np.array([q11DD, q21DD])
    qdDD=(A @ xiDD - B @ qDD - d)/l
    q12DD=qdDD[0]
    q22DD=qdDD[1]
    return q12DD, q22DD


def dynamic_model(q11, q12, q21, q22, q11D, q12D, q21D, q22D):
    M = np.zeros((2,2))
    c = 0
    A,B=compute_A_B(q11, q12, q21, q22,True)
    Z= np.array([[ZZ1R,0],[0,ZZ2R]])
    d=compute_d(q11, q12, q21, q22,q11D, q12D, q21D, q22D, True)
    J=np.linalg.inv(A)@B
    
    M=Z+mR*np.transpose(J)@J
    c=mR*np.transpose(J)@np.linalg.inv(A)@d
    return M, c

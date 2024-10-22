from lab_amoro.parallel_robot import *
from lab_amoro.plot_tools import *
from biglide_models import *  # Modify this to use the biglide
import sys

gamma1=-1
Kd=10
Kp=10

def polynomial(tf, steps, deriv):
	t=np.linspace(0,tf,steps)
	if deriv==0:
		S=10*(t/tf)**3-15*(t/tf)**4+6*(t/tf)**5
	elif deriv==1:
		S=30/tf*(t/tf)**2-60/tf*(t/tf)**3+30/tf*(t/tf)**4
	elif deriv==2:
		S=60/(tf*tf)*t/tf-180/(tf*tf)*(t/tf)**2+120/(tf*tf)*(t/tf)**3
	return S
	
def cartesian_trajectory(tf, steps, pi, pf):
	p=pi+polynomial(tf,steps,0)*(pf-pi)
	v=polynomial(tf,steps,1)*(pf-pi)
	a=polynomial(tf,steps,2)*(pf-pi)
	return p,v,a


def main(args=None):
	# Initialize and start the ROS2 robot interface
	rclpy.init(args=args)
	robot = Robot("biglide")  # Modify this to use the biglide
	start_robot(robot)	

	# Prepare plots
	app = QtGui.QApplication([])
	scope_joint1 = Scope("Joint 1", -0.5, 1.5)
	scope_joint2 = Scope("Joint 2", -1.5, 1.5)

	# Create the trajectory as arrays in Cartesian space (position, velocity, acceleration)
	tf=20000
	steps=int(tf/10)
	
	
	pi=np.array([[0],[0]])
	pf=np.array([[0.08],[5]])
	
	vi=np.array([[0],[0]])
	vf=np.array([[0],[0]])
	
	ai=np.array([[0],[0]])
	af=np.array([[0],[0]])
	
	
	
	q=np.zeros((2,steps))
	phi=np.zeros((2,steps))
	qD=np.zeros((2,steps))
	phiD=np.zeros((2,steps))
	qDD=np.zeros((2,steps))
	phiDD=np.zeros((2,steps))
	
	P,V,A=cartesian_trajectory(tf, steps, pi, pf)
	print(P)
	
	for t in range(steps):
		#print(P[0,t])
		q[0,t],q[1,t] = igm(P[0,t],P[1,t], gamma1, gamma1)
		phi[0,t],phi[1,t] = dgm_passive(q[0,t], q[1,t], gamma1)
		
		qD[0,t], qD[1,t] = ikm(q[0,t], phi[0,t], q[1,t], phi[1,t], V[0,t], V[1,t])
		phiD[0,t],phiD[1,t] = dkm_passive(q[0,t], phi[0,t], q[1,t], phi[1,t], qD[0,t], qD[1,t], V[0,t], V[1,t])
		
		qDD[0,t], qDD[1,t]= ikm2(q[0,t], phi[0,t], q[1,t], phi[1,t], qD[0,t], phiD[0,t], qD[1,t], phiD[1,t], A[0,t], A[1,t])

		
		

	# Create the trajectory as arrays in joint space using the inverse models (position, velocity, acceleration)

	index = 0

	# Controller
	try:
		robot.apply_efforts(0.0, 0.0)  # Required to start the simulation
		while True:
			if robot.data_updated():
				# Robot available data - This is the only data thet you can get from a real robot (joint encoders)
				q11 = robot.active_left_joint.position
				q21 = robot.active_right_joint.position
				q11D = robot.active_left_joint.velocity
				q21D = robot.active_right_joint.velocity
				

				# CTC controller
				print("qqq",q11,q21)
				q12,q22 = dgm_passive(q11, q21, gamma1)
				q12D,q22D = dkm_passive(q11, q12, q21, q22, q11D, q21D, 1, 1)
				
				M,c = dynamic_model(q11, q12, q21, q22, q11D, q12D, q21D, q22D)
				
				qa= np.array([q11,q21])
				qaD=np.array([q11D,q21D])
				
				alfa=qDD[:,index] + Kd*( qD[:,index] - q11D) + Kp*(q[:,index] - qa)
				tau=np.matmul(M,alfa)+c
				
				M0 = np.array([M[0,0],M[0,1]])
				M1 = np.array([M[1,0],M[1,1]])
				tau_left=tau[0]
				 
				#np.dot(M0,qDD[0,index] + Kd*( qD[0,index] - q11D) + Kp*(q[0,index] - q11)) + c[0]
				tau_right=tau[1] 
				#np.dot(M1,qDD[1,index] + Kd*( qD[1,index] - q21D) + Kp*(q[1,index] - q21)) + c[1]
				
				robot.apply_efforts(tau_left, tau_right)
				

				# Scope update
				time = robot.get_time()
				if time < 5.0:
					scope_joint1.update(time, 0.0, 0.0)
					scope_joint2.update(time, 0.0, 0.0)

				if index < steps-1:
					index += 1  # Next point in trajectory

	except KeyboardInterrupt:
		pass


if __name__ == "__main__":
	main(sys.argv)

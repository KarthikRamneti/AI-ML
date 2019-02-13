import cmath
import numpy as np
import sympy
import matplotlib.pyplot as plt
import subprocess
import shlex

def dir_vec(AB):
	return np.matmul(AB,dvec)

def norm_vec(AB):
	return np.matmul(omat,np.matmul(AB,dvec))

def line_intersect(AD,CF):
	n1 = norm_vec(AD)
	n2 = norm_vec(CF)
	N = np.vstack((n1,n2))
	p = np.zeros(2)
	p[0] = np.matmul(n1,AD[:,0])
	p[1] = np.matmul(n2,CF[:,0])
	return np.matmul(np.linalg.inv(N),p)

def altitude_side_intersection_pt(BC,A,B):
	n1 = dir_vec(BC)
	n2 = norm_vec(BC)
	N = np.vstack((n1,n2))
	p = np.zeros(2)
	p[0] = np.matmul(n1,A)
	p[1] = np.matmul(n2,B)
	return np.matmul(np.linalg.inv(N),p)

def Solve_Quadratic(A):
	a = A[0]
	b = A[1]
	c = A[2]
	d = (b**2) - (4*a*c)
	x = (-b-cmath.sqrt(d))/(2*a)
	y = (-b+cmath.sqrt(d))/(2*a)
	if d >= 0:
		return np.array([float(x),float(y)])
	else:
		return 0.01

Area = 28

x = sympy.Symbol('k')
m = sympy.Matrix([[x,-3*x,1],[5,x,1],[-1*x,2,1]])
Det = m.det()

Coeff_det = m.det().as_poly().coeffs()
Coeff_eqn1 = m.det().as_poly().coeffs()
Coeff_eqn1[-1] = Coeff_det[-1] - 2*Area
k1 = Solve_Quadratic(Coeff_eqn1)

Coeff_eqn2 = m.det().as_poly().coeffs()
Coeff_eqn2[-1] = Coeff_det[-1] + 2*Area
k2 = Solve_Quadratic(Coeff_eqn2)

k= np.hstack((k1,k2))
l = len(k)
for i in range(0,l):
	X = k[i]%1
	if X == 0 :
		k_value = k[i]
		#print k[i]
		
A = np.array([k_value,-3*k_value])
B = np.array([5,k_value])
C = np.array([-1*k_value,2])

AB = np.vstack((A,B)).T
BC = np.vstack((B,C)).T
CA = np.vstack((C,A)).T

dvec = np.array([-1,1])
omat = np.array([[0,1],[-1,0]])

P = altitude_side_intersection_pt(BC,A,B)
Q = altitude_side_intersection_pt(CA,B,C)
R = altitude_side_intersection_pt(AB,C,A)


AP = np.vstack((A,P)).T
BQ = np.vstack((B,Q)).T
CR = np.vstack((C,R)).T

H1 = line_intersect(AP,CR)
H2 = line_intersect(AP,BQ)
H3 = line_intersect(BQ,CR)

len = 10
lam_1 = np.linspace(0,1,len)

x_AB = np.zeros((2,len))
x_BC = np.zeros((2,len))
x_CA = np.zeros((2,len))
x_AP = np.zeros((2,len))
x_BQ = np.zeros((2,len))
x_CR = np.zeros((2,len))
for i in range(len):
	temp1 = A + lam_1[i]*(B-A)
	x_AB[:,i]= temp1
	temp2 = B + lam_1[i]*(C-B)
	x_BC[:,i]= temp2
	temp3 = C + lam_1[i]*(A-C)
	x_CA[:,i]= temp3
	temp4 = A + lam_1[i]*(P-A)
	x_AP[:,i]= temp4
	temp5 = B + lam_1[i]*(Q-B)
	x_BQ[:,i]= temp5
	temp6 = C + lam_1[i]*(R-C)
	x_CR[:,i]= temp6

plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_AP[0,:],x_AP[1,:],label='$AP$')
plt.plot(x_BQ[0,:],x_BQ[1,:],label='$BQ$')
plt.plot(x_CR[0,:],x_CR[1,:],label='$CR$')

plt.plot(A[0],A[1], 'o')
plt.text(A[0]*(1+0.05),A[1]*(1+0.01),'A')
plt.plot(B[0],B[1], 'o')
plt.text(B[0]*(1+0.03),B[1]*(1),'B')
plt.plot(C[0],C[1], 'o')
plt.text(C[0]*(1+0.1),C[1]*(1+0.03),'C')
plt.plot(P[0],P[1], 'o')
plt.text(P[0]*(1+0.1),P[1]*(1+0.03),'P')
plt.plot(Q[0],Q[1], 'o')
plt.text(Q[0]*(1+0.5),Q[1]*(1-0.1),'Q')
plt.plot(R[0],R[1], 'o')
plt.text(R[0]*(1+0.04),R[1]*(1),'R')
plt.plot(H1[0],H1[1], 'o')
plt.text(H1[0]*(1+0.15),H1[1]*(1-0.27),'H')

plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid()

plt.show()

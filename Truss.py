import sys
import numpy as np 
import math
from matplotlib import pyplot as plt
import copy

# Constants
E = 200*10**9 #[N/mm2]
A = 0.005 #[m2]
nDoF = 6 # Total number of degrees of freedom
restaindDoF = [1,2,5,6]

# Function to calculate the member global stiffness matrix
def calculateKg(E,A,L,theta):
    c = math.cos(theta)
    s = math.sin(theta) 
    


    K11 = (E*A/L)*np.array([[c**2,c*s],[c*s,s**2]])
    K12 = (E*A/L)*np.array([[-c**2,-c*s],[-c*s,-s**2]])
    K21 = (E*A/L)*np.array([[-c**2,-c*s],[-c*s,-s**2]])
    K21 = (E*A/L)*np.array([[c**2,c*s],[c*s,s**2]])

    return[K11,K12,K21,K21]


# 1. Global stiffness matrix for each member
[K11_12, K12_12, K21_12, K22_12] =  calculateKg(E,A,3,0) # Member A
[K11_23, K12_23, K21_23, K22_23] =  calculateKg(E,A,5,2.2143) # Member B


# 2. Primary stiffness matrix for the structure
k11 = K11_12
k12 = K12_12
k13 = np.zeros([2,2])

k21 = K21_12
k22 = K22_12 + K11_23
k23 = np.zeros([2,2])

k31 = np.zeros([2,2])
k32 = K12_12
k33 = K22_23

r1 = np.concatenate((k11,k12,k13),axis=1)
r2 = np.concatenate((k21,k22,k23),axis=1)
r3 = np.concatenate((k31,k32,k33),axis=1)
Kp = np.concatenate((r1,r2,r3),axis=0)

# 3. Extract the structure stiffness matrix, Ks
restrainedIndex = [x-1 for x in restaindDoF] #Index for each restrained DoF (lists starts at 0)

# Reduce to structure stiffness matrix by deleting thew rows in and columns for restrained DoF
Ks = np.delete(Kp,restrainedIndex, 0) # Delete rows
Ks = np.delete(Ks,restrainedIndex, 1) # Delete columns
Ks = np.matrix(Ks) # convert from Numpy array to a matrix

# 4. Solve unknown displacements
U2 = Ks.I*np.array([[0],[-150000]])
U_x2 = U2[0,0]
U_y2 = U2[1,0]


# 5. Determine the reaction forces
UG = np.array([[0,0,U_x2, U_y2,0,0]])    # Global displ vector
FG = np.matmul(Kp,UG.T)                  # Force vector. Use the global primary stiffness matrix

# 6. Solve for member forces 
# Write a function that calculates member forces based on nodal displacements
def calculateForce(E,A,L,theta, Ua,Ub):

    #Transformation matrix
    c = math.cos(theta)
    s = math.sin(theta)
    T = np.array([[c,s,0,0],[0,0,c,s]]) # Transformation matrix - Global to local
    disp = np.array([[Ua[0], Ua[1], Ub[0], Ub[1]]]).T # Global displacements
    disp_local = np.matmul(T,disp) # Local displacements
    F_axial = (E*A/L)*(disp_local[1] - disp_local[0]).item()
    return F_axial

F12 = calculateForce(E,A,3,0,[0,0],[U_x2,U_y2])
F23 = calculateForce(E,A,5,2.2143,[U_x2,U_y2],[0,0])

print(F23)

# 7. Plotting
# Scale displacement
xFac = 100


fig = plt.figure()
axes = fig.add_axes([0,0,1,1])
fig.gca().set_aspect('equal', adjustable ='box')

## Structure
axes.plot([0,3],[0,0],'b')
axes.plot([0,3],[4,0],'b')

# Deformed structure 
axes.plot([0,3+U_x2*xFac],[0,0+U_y2*xFac],'--r')
axes.plot([0,3+U_x2*xFac],[4,0+U_y2*xFac],'--r')

# Plot points
axes.plot([0],[0],'bo')
axes.plot([3],[0],'bo')
axes.plot([0],[4],'bo')

# Plot deflection label
plt.text(3+U_x2*xFac,0+U_y2*xFac, 'Ux = {one}, Uy = {two}'.format(one = np.round(U_x2,4), two = np.round(U_y2,4)))

axes.grid()
axes.set_xlabel('Distance [m]')
axes.set_ylabel('Distance [m]')
axes.set_title('Deflected shape')


plt.show()
print(U_x2)

print(fig)

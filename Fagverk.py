import sys
import numpy as np 
import math
from matplotlib import pyplot as plt
import copy

# INPUT DATA
# -------------------
# Constants
E = 200*10**9  #[N/m2]
A = 0.005 # [m2]
xFac = 100 # Scale factor for plotted displacement

# Nodal coordinates [x,y], in ascending node order
nodes = np.array([[0,0],
                [4,0],
                [8,0],
                [12,0],
                [10,3],
                [6,3],
                [2,3]])

# Members [node_i, node_j]
members = np.array([[1,2],
                   [2,3],
                   [3,4],
                   [4,5],
                   [5,6],
                   [6,7],
                   [7,1],
                   [7,2],
                   [2,6],
                   [6,3],
                   [3,5]])

# Supports
restrainedDoF = [1,2,10] # The degrees of freedom restrained by supports

# Loading
#forceVector = np.array([[0,-200000,0,0,0,0,0,0,0,0,0,0]]).T # Vector of externally applied loads
forceVector = np.array([np.zeros(14)]).T
#forceVector[3] = -500000
forceVector[5] = -500000
#forceVector[7] = -100000
# END OF DATA ENTRY
# -------------------



# Calculate member orientation and length
#----------------------------------------
# Define a function to calculate the member orientation and length
def memberOrientation(memberNo):
    memberIndex = memberNo-1 #Index identifying member in array of members
    node_i = members[memberIndex][0] # First node of the member
    node_j = members[memberIndex][1] # Second node of the member

    xi = nodes[node_i-1][0] # X-coord of node i
    yi = nodes[node_i-1][1] # Y-coord of node i

    xj = nodes[node_j-1][0] # X-coord of node j
    yj = nodes[node_j-1][1] # Y-coord of node j

    # Angle of member with respect to horizontal axis

    dx = xj-xi #x-component of member
    dy = yj-yi #y-component of member

    mag = math.sqrt(dx**2+dy**2) # Magnitude of vector (length of member)
    memberVector = np.array([dx,dy]) # Member represented as a vector

    # Need to capture quadrant first, then identify appropriate reference axis and offet angle
    if (dx>0 and dy==0):                                    #X-axis pos
        theta = 0
    elif (dx==0 and 0<dy):                                  #Y-axis pos
        theta = math.pi/2
    elif (dx<0 and dy==0):                                  #X-axis neg
        theta = math.pi
    elif (dx==0 and dy<0):                                  #Y-axis neg
        theta = 3*math.pi/2

    elif (0<dx and 0<dy):            
        #0<theta<90                   
        refVector = np.array([1,0])                         #Vector describing the positive x-axis
        theta = math.acos(refVector.dot(memberVector)/(mag))  #Standard formula for angle between 2 vectors

    elif (dx<0 and dy>0):  
        #90<theta<180                                 
        refVector = np.array([0,1])                         #Vector describing the positive y-axis
        theta = math.acos(refVector.dot(memberVector)/(mag)) + (math.pi/2)

    elif (dx<0 and dy<0):         
        #180<theta<270                         
        refVector = np.array([-1,0])                         #Vector describing the positive y-axis
        theta = math.acos(refVector.dot(memberVector)/(mag)) + (math.pi)

    elif (0<dx and dy<0):            
        #270<theta<360                       
        refVector = np.array([0,-1])                         #Vector describing the positive y-axis
        theta = math.acos(refVector.dot(memberVector)/(mag)) + (3*math.pi/2)

    return [theta, mag]

# Calculate orientation and length for each member and store
orientations = np.array([]) # Initialise an array to hold orientations
lengths  = np.array([]) # Initialise an array to hold member lengths

for n, mbr in enumerate(members):
    [angle, length] = memberOrientation(n+1) # Member 1, not index 0, because first n = 0
    orientations = np.append(orientations,angle)
    lengths = np.append(lengths,length)

# Define a function to calculate member global stiffness matrix
#--------------------------------------------------------------
def calculateKg(memberNo):
    theta = orientations[memberNo-1]
    mag = lengths[memberNo-1]

    c = math.cos(theta)
    s = math.sin(theta) 


    K11 = (E*A/mag)*np.array([[c**2,c*s],[c*s,s**2]])
    K12 = (E*A/mag)*np.array([[-c**2,-c*s],[-c*s,-s**2]])
    K21 = (E*A/mag)*np.array([[-c**2,-c*s],[-c*s,-s**2]])
    K22 = (E*A/mag)*np.array([[c**2,c*s],[c*s,s**2]])

    return[K11, K12, K21, K22]

# Build primart stiffness matrix, Kp
#-----------------------------------
nDoF = np.amax(members)*2 # Total numbers of degrees of freedom in the problem
Kp = np.zeros([nDoF,nDoF]) # Initialising the primary stiffness matrix

for m, mbr in enumerate(members):
    [K11, K12, K21, K22] = calculateKg(m+1)

    node_i = mbr[0] # Node number for node i of this member
    node_j = mbr[1] # Node number for node j of this member

    ia = 2*node_i-2
    ib = 2*node_i-1
    ja = 2*node_j-2
    jb = 2*node_j-1

    Kp[ia:ib+1,ia:ib+1] = Kp[ia:ib+1,ia:ib+1] + K11
    Kp[ia:ib+1,ja:jb+1] = Kp[ia:ib+1,ja:jb+1] + K12
    Kp[ja:jb+1,ia:ib+1] = Kp[ja:jb+1,ia:ib+1] + K21
    Kp[ja:jb+1,ja:jb+1] = Kp[ja:jb+1,ja:jb+1] + K22

# Extract the structure stiffness matrix, Ks
#---------------------------------------
restrainedIndex = [x-1 for x in restrainedDoF] #Index for each restrained DoF (lists starts at 0)

# Reduce to structure stiffness matrix by deleting the rows and columns for restrained DoF
Ks = np.delete(Kp,restrainedIndex, 0) # Delete rows
Ks = np.delete(Ks,restrainedIndex, 1) # Delete columns
Ks = np.matrix(Ks) # convert from Numpy array to a matrix


# Solve unknown displacements
# -------------------------------
forcevectorRed = copy.copy(forceVector)
forcevectorRed = np.delete(forcevectorRed, restrainedIndex, 0)
U = Ks.I*forcevectorRed


# Solve for reactions
# -------------------------------
#Construct the global displacement vector
UG = np.zeros(nDoF)

c = 0
for i in np.arange(nDoF):
    if i in restrainedIndex:
        UG[i] = 0
    else:
        UG[i] = U[c]
        c=c+1

UG = np.array([UG]).T
FG = np.matmul(Kp,UG)

# Solve for member forces
#--------------------------
mbrForces = np.array([]) # Initialize an array to hold member forces

for n, mbr in enumerate(members):
    theta = orientations[n]
    mag = lengths[n]

    node_i = mbr[0]
    node_j = mbr[1]

    ia = 2*node_i-2
    ib = 2*node_i-1
    ja = 2*node_j-2
    jb = 2*node_j-1

    # Transformation matrix
    c = math.cos(theta)
    s = math.sin(theta) 
    T = np.array([[c,s,0,0],[0,0,c,s]])

    disp = np.array([[ UG[ia], UG[ib], UG[ja], UG[jb] ]]).T
    disp_local = np.matmul(T,disp)[0]
   
    F_axial = (A*E/mag)*(disp_local[1]-disp_local[0]) # Axial load
    mbrForces = np.append(mbrForces,F_axial)

    
# Plotting
fig = plt.figure()
axes = fig.add_axes([0,0,1,1])
fig.gca().set_aspect('equal', adjustable ='box')

# Plot members
for mbr in members:
    node_i = mbr[0] # Node number for node i of this member
    node_j = mbr[1] # Node number for node j of this member

    ix = nodes[node_i-1,0] # x-coord of node i of this member
    iy = nodes[node_i-1,1] # y-coord of node i of this member
    jx = nodes[node_j-1,0] # x-coord of node j of this member
    jy = nodes[node_j-1,1] # y-coord of node j of this member

    # Index for DoF for this member
    ia = 2*node_i-2 # Horisozontal DoF at node i for this member
    ib = 2*node_i-1 # Vertical DoF at node i for this member
    ja = 2*node_j-2 # Horisozontal DoF at node j for this member
    jb = 2*node_j-1 # Vertical DoF at node j for this member

    axes.plot([ix,jx],[iy,jy],'b')
    axes.plot([ix + UG[ia,0]*xFac, jx + UG[ja,0]*xFac], [iy + UG[ib,0]*xFac,jy + UG[jb,0]*xFac],'--r')

# Plot nodes 
for node in nodes:
    axes.plot([node[0]],[node[1]],'bo')

axes.set_xlabel('Lengde [m]')
axes.set_ylabel('Lengde [m]')
axes.set_title('Nedbøyning')
axes.legend('upper right',['udeformert konstruksjon'])


plt.show()
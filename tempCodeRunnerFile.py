import os
import math
import gmsh
import platform
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eig

matplotlib.use('Agg')  # Use a non-interactive backend

################################################################################################
#                                    Geometry + Main Driver                                    #
################################################################################################

def main():
    ## Define geometry
    L = 125  # Length of each beam
    thetas = np.array([30])  # degrees
    elements_per_line = 10   # Elements per beam

    ## Material + parameters
    E = 210e3
    nu = 0.3
    beta = 1/15
    b = beta * L
    h = beta * L
    A = b * h
    I = (b * h**3)/12
    ngp = 2
    rho = 7850
    omega0 = 1e6*((np.pi/L)**2)*np.sqrt((E*I)/(rho*A))
    omeg_mult = 10

    mesh_type = "curved"  # <=== Change to "curved" to use curved geometry

    for theta_deg in thetas:
        theta = math.radians(theta_deg)

        if mesh_type == "straight":
            create_mesh(L, theta, elements_per_line)
        elif mesh_type == "curved":
            create_mesh_curved_beam(L, theta, elements_per_line)
        else:
            raise ValueError("Invalid mesh type!")

        conn_mat = np.loadtxt('Data Files/conn_mat.dat', dtype=int) - 1
        coordinates = np.loadtxt('Data Files/coordinates.dat')
        const_param = np.array([E, nu, A, I, ngp, rho, omega0, L, theta, omeg_mult])
        analysis(const_param, conn_mat, coordinates, theta, theta_deg)

# ------------------------------------------------------------------------------------------------
# Straight Beam Mesh Function
# ------------------------------------------------------------------------------------------------
def create_mesh(L, theta, elements_per_line):
    if not os.path.exists('Data Files'):
        os.makedirs("Data Files")

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("unit_cell")

    p1 = (0.0, 0.0, 0.0)
    p2 = (0.0, L, 0.0)
    p3 = (L * math.cos(theta), L * (1 + math.sin(theta)), 0.0)
    p4 = (-L * math.cos(theta), L * (1 + math.sin(theta)), 0.0)

    coords = []
    conn = []
    coord_to_tag = {}

    def add_point(pt):
        if pt not in coord_to_tag:
            tag = gmsh.model.geo.addPoint(*pt, 0)
            coord_to_tag[pt] = tag
            coords.append((pt[0], pt[1]))
        return coord_to_tag[pt]

    def linspace_points(pa, pb, n):
        return [
            (
                pa[0] + i * (pb[0] - pa[0]) / n,
                pa[1] + i * (pb[1] - pa[1]) / n,
                0.0
            ) for i in range(n)
        ] + [pb]

    beam1 = linspace_points(p1, p2, elements_per_line)
    beam2 = linspace_points(p2, p4, elements_per_line)[1:]
    beam3 = linspace_points(p2, p3, elements_per_line)[1:]

    prev = add_point(beam1[0])
    for pt in beam1[1:]:
        curr = add_point(pt)
        gmsh.model.geo.addLine(prev, curr)
        conn.append((prev, curr))
        prev = curr

    center_id = prev
    beam2_tags = []
    beam3_tags = []
    b2_idx, b3_idx = 0, 0
    len_b2, len_b3 = len(beam2), len(beam3)
    turn = 0

    while b2_idx < len_b2 or b3_idx < len_b3:
        if turn == 0 and b2_idx < len_b2:
            tag = add_point(beam2[b2_idx])
            beam2_tags.append(tag)
            b2_idx += 1
        elif turn == 1 and b3_idx < len_b3:
            tag = add_point(beam3[b3_idx])
            beam3_tags.append(tag)
            b3_idx += 1
        turn = 1 - turn

    conn.append((center_id, beam2_tags[0]))
    gmsh.model.geo.addLine(center_id, beam2_tags[0])
    for i in range(len(beam2_tags) - 1):
        conn.append((beam2_tags[i], beam2_tags[i + 1]))
        gmsh.model.geo.addLine(beam2_tags[i], beam2_tags[i + 1])

    conn.append((center_id, beam3_tags[0]))
    gmsh.model.geo.addLine(center_id, beam3_tags[0])
    for i in range(len(beam3_tags) - 1):
        conn.append((beam3_tags[i], beam3_tags[i + 1]))
        gmsh.model.geo.addLine(beam3_tags[i], beam3_tags[i + 1])

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(1)

    tag_list = list(coord_to_tag.values())
    tag_to_index = {tag: i + 1 for i, tag in enumerate(tag_list)}

    with open('Data Files/coordinates.dat', "w") as f:
        for tag in tag_list:
            x, y, _ = gmsh.model.getValue(0, tag, [])
            f.write(f"{x:.6f}\t{y:.6f}\n")

    with open('Data Files/conn_mat.dat', "w") as f:
        for n1, n2 in conn:
            f.write(f"{tag_to_index[n1]}\t{tag_to_index[n2]}\n")

    gmsh.fltk.run()
    gmsh.finalize()

# ------------------------------------------------------------------------------------------------
# Curved Beam Mesh Function
# ------------------------------------------------------------------------------------------------
def create_mesh_curved_beam(L, theta, elements_per_line):
    if not os.path.exists('Data Files'):
        os.makedirs("Data Files")

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("curved_Y_beam")

    # Beam parameters
    radius = L / theta
    dtheta = theta / elements_per_line

    coords = []
    conn = []
    coord_to_tag = {}

    def add_point(pt):
        if pt not in coord_to_tag:
            tag = gmsh.model.geo.addPoint(*pt, 0)
            coord_to_tag[pt] = tag
            coords.append((pt[0], pt[1]))
        return coord_to_tag[pt]

    # Vertical straight beam (p1 → p2)
    p1 = (0, 0, 0)
    p2 = (0, L, 0)
    vert_points = [(0, i * L / elements_per_line, 0) for i in range(elements_per_line + 1)]
    vert_tags = []

    prev = add_point(vert_points[0])
    for pt in vert_points[1:]:
        curr = add_point(pt)
        gmsh.model.geo.addLine(prev, curr)
        conn.append((prev, curr))
        vert_tags.append(curr)
        prev = curr

    top_center = vert_tags[-1]

    # Curved left arm (counter-clockwise)
    left_points = [
        (
            -radius * math.sin(i * dtheta),
            L + radius * (1 - math.cos(i * dtheta)),
            0.0
        )
        for i in range(1, elements_per_line + 1)
    ]

    prev = top_center
    for pt in left_points:
        curr = add_point(pt)
        gmsh.model.geo.addLine(prev, curr)
        conn.append((prev, curr))
        prev = curr

    # Curved right arm (clockwise)
    right_points = [
        (
            radius * math.sin(i * dtheta),
            L + radius * (1 - math.cos(i * dtheta)),
            0.0
        )
        for i in range(1, elements_per_line + 1)
    ]

    prev = top_center
    for pt in right_points:
        curr = add_point(pt)
        gmsh.model.geo.addLine(prev, curr)
        conn.append((prev, curr))
        prev = curr

    # Finalize
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(1)

    tag_list = list(coord_to_tag.values())
    tag_to_index = {tag: i + 1 for i, tag in enumerate(tag_list)}

    # Write coordinates.dat
    with open('Data Files/coordinates.dat', 'w') as f:
        for tag in tag_list:
            x, y, _ = gmsh.model.getValue(0, tag, [])
            f.write(f"{x:.6f}\t{y:.6f}\n")

    # Write conn_mat.dat
    with open('Data Files/conn_mat.dat', 'w') as f:
        for n1, n2 in conn:
            f.write(f"{tag_to_index[n1]}\t{tag_to_index[n2]}\n")

    gmsh.fltk.run()
    gmsh.finalize()


    
# Function: Perform analysis
def analysis(const_param, conn_mat, coordinates, theta, theta_d):
    # Unpack all constant parameters
    E, nu, A, I, ngp, rho, omega0, L, theta, omeg_mult = const_param
    ngp, rho, omeg_mult = int(ngp), int(rho), int(omeg_mult)
    
    nel, nnode = len(conn_mat), len(coordinates)                         # Nunmber of elements and nodes
    print(f'''_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _\n'''+
          f'''\nThe given finite element setup has {nnode} nodes and {nel} elements'''+
          '''\n_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _\n''')
        
    ## Discretisation
    x, y = coordinates[:, 0].reshape(-1, 1), coordinates[:, 1].reshape(-1, 1)        # Column vectors --> x and y coordinates of nodes
    dof = 3                                     # Degrees of freedom per node 
    ndof = dof*nnode                            # Total degrees of freedom in the system
    
    ## Solution
    Axx, Dxx = E*A, E*I                         # Extension and bending stiffness coefficients
    gp, gw = gaussLegendreQuadrature(ngp, 0, 1)         # Get Gauss-Legendre quadrature points and weights (linear stiffness terms)
    
    print('Forming mass matrix')
    M = mass_mat(rho*A, 0, rho*I, ndof, conn_mat, x, y, ngp, nel, gp, gw)
    print('Mass matrix formation complete\n_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _\n')
    print('Forming stiffness matrix')
    K = stiffness_mat(ngp, ndof, nel, conn_mat, x, y, gp, gw, Axx, Dxx)
    print('Stiffness matrix formation complete\n_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _\n')
    
    print('Performing dispersion analysis')
    dispersion(K, M, omega0, 0, 3, nnode, theta, theta_d, omeg_mult)
    print('Dispersion analysis complete\n_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _\n')
    
# Function: Gauss-Legendre Quadrature points and weights
def gaussLegendreQuadrature(n, a, b):
    beta = 0.5 / np.sqrt(1 - (2*np.array(list(range(1, n)), dtype=float))**(-2))        # 3-term recurrence coefficients
    J = np.diag(beta, 1) + np.diag(beta, -1)                                            # Jacobi matrix
    zeta, V = np.linalg.eigh(J)                                                         # Eigenvalue decomposition
    zeta, i = np.sort(zeta), np.argsort(zeta)                                           # Sorted zeta and indices for sorted zeta
    w = 2 * V[0, i]**2                                                                  # Quadrature weights
    
    # Scale points from [-1, 1] to [a, b]
    zeta = (b-a)/2 * zeta + (a+b)/2
    w = (b-a)/2*w
    
    return zeta, w

# Function: Mass Matrix
def mass_mat(I0, I1, I2, ndof, conn_mat, x, y, ngp, nel, gp, gw):
    # Initialization for lobal assembly
    M = np.zeros((ndof, ndof))  # Initialize global mass matrix
        
    for globitr in range(nel):  # Loop: global assembly --> start
        M11, M12, M21, M22 = np.zeros((2, 2)), np.zeros((2, 4)), np.zeros((4, 2)), np.zeros((4, 4))  # Initialize mass submatrices M11, M12, M21 and M22
        
        curr_dof = element_dof(conn_mat[globitr, :], 3)      # Extracts degrees of freedom associated with the current element        
        detJ = np.sqrt( (y[conn_mat[globitr, 1], 0] - y[conn_mat[globitr, 0], 0])**2 + (x[conn_mat[globitr, 1], 0] - x[conn_mat[globitr, 0], 0])**2 )   # Determinant of Jacobian
        
        for glitr in range(ngp):    # Loop: numerical integration --> start
            zeta, w = gp[glitr], gw[glitr]        # Current Gauss point and Gauss weight
            
            # Relevenat shape functions and/or derivatives
            Nu, _ = axial_shape_func(zeta, detJ)
            Nw, dNwdx, _, _ = transverse_shape_func(zeta, detJ)
            Nu = np.array(Nu).reshape(2, 1)
            Nw = np.array(Nw).reshape(1, 4)
            dNwdx = np.array(dNwdx).reshape(1, 4)
            
            # Mass submatrix M11 (2X2)
            M11 += w*detJ*I0*(Nu @ Nu.T)
            
            # Mass submatrix M12 (2X4)
            M12 -= w*detJ*I1*(Nu @ dNwdx)
            
            # Mass submatrix M21 (4X2)
            M21 -= w*detJ*I1*(dNwdx.T @ Nu.T)
            
            # Mass submatrix M22 (4X4)
            M22 += w*detJ*( I0* (Nw.T @ Nw) + I2*(dNwdx.T @ dNwdx) )
                                
        # Loop: numerical integration --> end
        
        Me = np.block([[M11, M12],
                      [M21, M22]])
        
        # Matrix rearrangement
        tmp = Me[:, [1]]
        Me[:, [1]], Me[:, [2]] = Me[:, [2]], Me[:, [3]]
        Me[:, [3]] = tmp
        tmp = Me[[1], :]
        Me[[1], :], Me[[2], :] = Me[[2], :], Me[[3], :]
        Me[[3], :] = tmp

        Me = coord_trans(Me, x[conn_mat[globitr, :]], y[conn_mat[globitr, :]], 'l2g')   # Coordinate transformation        
        
        M[np.ix_(curr_dof, curr_dof)] += Me      # Assembly to current (working) degrees of freedom
        
    # Loop: global assembly --> end
    M = M.astype('float32')
    
    return M

# Function: Gets degrees of freedom associated with current element
def element_dof(nodes, ndof):
    edof = np.zeros(ndof) if isinstance(nodes, (int, float, complex)) else np.zeros(ndof * nodes.size, dtype=int)
    n = 1 if isinstance(nodes, (int, float, complex)) else nodes.size
    
    for i in range(n):
        base_index = 3 * (nodes[i]+1)
        edof[3*i:3*i+3] = np.array([base_index - 2, base_index - 1, base_index]) - 1

    return edof

# Function: Linear shape functions (axial dofs) and relevant derivatives
def axial_shape_func(zeta, l):
    Nu = np.array([1-zeta, zeta])       # Shape functions
    dNudx = (1/l)*np.array([-1, 1])     # First derivatives
    
    return Nu, dNudx

# Function: Hermite cubic polynomial shape functions (transverse and rotational dofs) and relevant derivatives
def transverse_shape_func(zeta, l):
    Nw = np.array([1 - 3*zeta**2 + 2*zeta**3,
                   -l*zeta*(1 - zeta)**2,
                   3*zeta**2 - 2*zeta**3,
                   -l*zeta*(zeta**2 - zeta)])               # Shape functions
    dNwdx = (1/l)*np.array([-6*zeta + 6*zeta**2,
                            -l*(1 + 3*zeta**2 - 4*zeta),
                            6*zeta - 6*zeta**2,
                            -l*(3*zeta**2 - 2*zeta)])       # First derivatives
    d2Nwdx2 = (1/l**2)*np.array([12*zeta - 6,
                                 -l*(6*zeta - 4),
                                 6 - 12*zeta,
                                 -l*(6*zeta - 2)])          # Second derivatives
    d3Nwdx3 = (1/l**3)*np.array([12, -6*l, 12, -6*l])       # Third derivatives
    
    return Nw, dNwdx, d2Nwdx2, d3Nwdx3

# Function: Coordinate transformation of matrices and vectors
def coord_trans(A, x, y, tr_case):
    c, s = trig_dat(x, y)
    Q0 = np.array([[c, s, 0],
                   [-s, c, 0],
                   [0, 0, 1]])
    Q = np.block( [ [Q0, np.zeros((3, 3))],
                    [np.zeros((3, 3)), Q0]] )       # Coordinate transformation matrix
    
    s_info = A.shape
    
    if tr_case == 'l2g':                            # Transformation: Local to Global
        if s_info[1]!=1:                            # Check if A has only 1 column (vector) or not (matrix)
            Tr = np.dot(np.dot((Q.T), A), Q)        # Coordinate transformation for stiffness and mass matrices
        else:
            Tr = np.dot(Q.T, A)                       # Coordinate transformation for force and displacement vectors
    elif tr_case == 'g2l':                          # Transformation: Global to Local
        if s_info[1]!=1:                            # Check if A has only 1 column (vector) or not (matrix)
            Tr = np.dot( np.dot(Q, A), Q.T)        # Coordinate transformation for stiffness and mass matrices
        else:
            Tr = np.dot(Q, A)                     # Coordinate transformation for force and displacement vectors
    else:
        raise ValueError('Invalid Coordinate Transformation!!')
        
    return Tr

# Function: Angle
def trig_dat(x, y):
    x1, x2, y1, y2 = x[0, 0], x[1, 0], y[0, 0], y[1, 0]
    alpha = np.arctan2(y2-y1, x2-x1)
    if alpha < 0:
        alpha += 2*np.pi
    [c, s] = [np.cos(alpha), np.sin(alpha)]
    return c, s

# Function: Stiffness Matrix, Tangent Stiffness MAtrix, Force Vector and Residual Force Vector
def stiffness_mat(ngp_curr, ndof, nel, conn_mat, x, y, gp, gw, Axx, Dxx):
    K = np.zeros((ndof, ndof))                      # Initialize global stiffness matrix
    
    for globitr in range(nel):                      # Loop: global assembly --> start
        K11, K12, K21, K22 = np.zeros((2, 2)), np.zeros((2, 4)), np.zeros((4, 2)), np.zeros((4, 4))  # Initialize stiffness matrix submatrices K11, K12, K21 and K22
        
        curr_dof = element_dof(conn_mat[globitr, :], 3)      # Extracts degrees of freedom associated with the current element
        
        detJ = np.sqrt( (y[conn_mat[globitr, 1], 0] - y[conn_mat[globitr, 0], 0])**2 + (x[conn_mat[globitr, 1], 0] - x[conn_mat[globitr, 0], 0])**2 )   # Determinant of Jacobian
         
        ## Linear stiffness coefficients and force vector
        for glitr in range(ngp_curr):            # Loop: numerical integration --> start
            zeta, w = gp[glitr], gw[glitr]                # Current Gauss point and Gauss weight
            
            # Relevant shape functions and/or derivatives
            Nu, dNudx = axial_shape_func(zeta, detJ)
            Nw, dNwdx, d2Nwdx2, _ = transverse_shape_func(zeta, detJ)
            
            # Reshape hape function vectors for outer products
            Nu = np.array(Nu).reshape(1, 2)
            dNudx = np.array(dNudx).reshape(1, 2)
            Nw = np.array(Nw).reshape(1, 4)
            dNwdx = np.array(dNwdx).reshape(1, 4)
            d2Nwdx2 = np.array(d2Nwdx2).reshape(1, 4)
            
            # Stiffness submatrix K11 (2X2)
            K11 += w*detJ*Axx* (dNudx.T @ dNudx)
            
            # Stiffness submatrix K22 (linear part) (4X4)
            K22 += w*detJ*Dxx* (d2Nwdx2.T @ d2Nwdx2)
            
        # Loop: numerical integration --> end
                
        Ke = np.block([[K11, K12],
                       [K21, K22]])
        
        
        # Matrix rearrangement:  Stiffness Matrix
        tmp = Ke[:, [1]]
        Ke[:, [1]], Ke[:, [2]] = Ke[:, [2]], Ke[:, [3]]
        Ke[:, [3]] = tmp
        tmp = Ke[[1], :]
        Ke[[1], :], Ke[[2], :] = Ke[[2], :], Ke[[3], :]
        Ke[[3], :] = tmp
        
        Ke = coord_trans(Ke, x[conn_mat[globitr, :]], y[conn_mat[globitr, :]], 'l2g')       # Coordinate transformation        
        
        ## Assembly into the global system
        # Assemble global stiffness matrix
        K[np.ix_(curr_dof, curr_dof)] += Ke

    # Loop: global assembly --> end
    
    K = K.astype('float32')
    
    return K

# Function: Dispersion analysis of unit cell    
def dispersion(K, M, omega0, case, ndof, nnode, theta, theta_d, omeg_mult):
    exit_nodes = 2
    M = M.astype('float64')*1e-12
    # Function: Identify the irreducible Brillouin zone
    def IBZ(theta):
        phi = np.pi/4 + theta/2        
        # Re-entrant
        if theta >= -np.pi/6 and theta < 0:
            Orecip = 2*np.pi*np.array([[0], [0]])
            Arecip = 2*np.pi*np.array([[1 / 2], [-1 / 2]])
            Brecip = 2*np.pi*np.array([[1 - 1 / (4 * (np.cos(phi))**2)], [-1 / (4 * (np.cos(phi))**2)]])
            Crecip = 2*np.pi*np.array([[1/(4*(np.cos(phi))**2)], [1/(4*(np.cos(phi))**2)]])
        # Honeycomb
        elif theta >= 0 and theta < np.pi/2:
            Orecip = 2*np.pi*np.array([[0], [0]])
            Arecip = 2*np.pi*np.array([[1 / (4*(np.sin(phi))**2)], [-1 / (4*(np.sin(phi))**2)]])
            Brecip = 2*np.pi*np.array([[1 - 1 / (4 * (np.sin(phi))**2)], [-1 / (4 * (np.sin(phi))**2)]])
            Crecip = 2*np.pi*np.array([[1/2], [1/2]])
        
        return Orecip.flatten(), Arecip.flatten(), Brecip.flatten(), Crecip.flatten()
    
    O, A, B, C = IBZ(theta)
    
    def linspace_segment(x, y, num):
        kx = np.linspace(x[0], x[1], num)
        ky = y[0] + ((y[1]-y[0])/(x[1]-x[0]))*(kx-x[0])
        return kx, ky

    # Create segments for Bloch wave vector k
    k_OAx, k_OAy = linspace_segment(np.array([O[0], A[0]]), np.array([O[1], A[1]]), 50)
    k_ABx, k_ABy = linspace_segment(np.array([A[0], B[0]]), np.array([A[1], B[1]]), 50)
    k_BCx, k_BCy = linspace_segment(np.array([B[0], C[0]]), np.array([B[1], C[1]]), 50)
    k_COx, k_COy = linspace_segment(np.array([C[0], O[0]]), np.array([C[1], O[1]]), 50)
    
    k_x = np.concatenate([k_OAx[:len(k_OAx)], k_ABx[:len(k_ABx)], k_BCx[:len(k_BCx)], k_COx])
    k_y = np.concatenate([k_OAy[:len(k_OAy)], k_ABy[:len(k_ABy)], k_BCy[:len(k_BCy)], k_COy])

    k = np.vstack((k_x, k_y)).astype('float64')         # Assemble Bloch wave vectors along the IBZ boundary
    
    freq_dimless = np.zeros((ndof * (nnode - exit_nodes), k.shape[1]))
    
    # Initialize translation matrix
    Tmat = np.zeros((ndof * nnode, ndof * (nnode - exit_nodes)), dtype=complex)     # Initialize as complex type matrix to handle complex values
    Tmat[:ndof * (nnode - exit_nodes), :ndof * (nnode - exit_nodes)] = np.eye(ndof * (nnode - exit_nodes))
    
    eye_dof = np.eye(ndof)
    OABCO = []

    for itr in range(k.shape[1]):
        k_curr = k[:, itr]
        k1, k2 = 1j*k_curr[0], 1j*k_curr[1]
        
        OABCO.append(np.sqrt(k[0, itr]**2 + k[1, itr]**2))
            
        Tmat[ndof * (nnode - exit_nodes):ndof * (nnode - exit_nodes) + ndof, :3] = eye_dof * np.exp(k2)
        Tmat[ndof * (nnode - exit_nodes) + ndof:ndof * (nnode - exit_nodes) + 2 * ndof, :3] = eye_dof * np.exp(k1)
            
        K_tr, M_tr = (Tmat.conj().T @ K) @ Tmat, (Tmat.conj().T @ M) @ Tmat
            
        eigenvalues = eig(K_tr, M_tr, right=False)
        eigenvalues_real = eigenvalues.real
        eigenvalues_real[eigenvalues_real < 0] = 0
            
        freq = np.sqrt(eigenvalues_real).real
        freq_dimless[:, itr] = np.sort(freq/omega0)
    
    real_freq = np.real(freq_dimless[:omeg_mult, :k.shape[1]])
    size_freq = real_freq.shape[0]

    max_freq, min_freq = np.zeros(size_freq), np.zeros(size_freq)

    # Find max and min frequency for each row (frequency mode)
    for freq_itr in range(size_freq):
        max_freq[freq_itr], min_freq[freq_itr] = np.max(real_freq[freq_itr, :]), np.min(real_freq[freq_itr, :])

    max_freq, min_freq = np.delete(max_freq, size_freq-1), np.delete(min_freq, 0)
 
    diff_freq = min_freq - max_freq

    mid_point_freq = 0.5 * (min_freq + max_freq)

    negative_indices = np.where(diff_freq < 0)[0]
    mid_point_freq, max_freq, min_freq = np.delete(mid_point_freq, negative_indices), np.delete(max_freq, negative_indices), np.delete(min_freq, negative_indices)

    bound_len = len(OABCO)

    plt.plot(np.arange(0, bound_len), real_freq[:omeg_mult, :k.shape[1]].T, '-b', linewidth=2)
    freq_data = real_freq[:omeg_mult, :k.shape[1]]
    np.savetxt('wave_vec.dat', np.arange(0, bound_len).T, fmt = "%f")
    np.savetxt(f'freq_data_{int(theta_d)}.dat', freq_data, fmt = "%f")
    
    indices = [0, 
               len(k_OAx) - 1, 
               len(k_OAx) + len(k_ABx) - 1, 
               len(k_OAx) + len(k_ABx) + len(k_BCx) - 1, 
               len(k_OAx) + len(k_ABx) + len(k_BCx) + len(k_COx) - 1]
    
    labels = ['Γ', "M'", 'X', 'M', 'Γ']
    
    plt.xlim([0, k.shape[1]])
    plt.ylim([0, omeg_mult])
    plt.xticks(ticks=indices, labels=labels)
    plt.yticks(np.linspace(0, omeg_mult, 6))
    plt.xlabel('Bloch wave vector \u03BA', fontsize=24, color='black', fontdict={'family': 'Times New Roman'})
    plt.ylabel('Normalized Frequency \u03A9', fontsize=24, color='black', fontdict={'family': 'Times New Roman'})

    plt.gca().tick_params(labelsize = 24)
    plt.gca().set_facecolor('white')
    
    output_dir = os.path.join('Figures', 'Dispersion')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    filename = (os.path.join(output_dir, 'LinD.pdf') if case == 0 else
                os.path.join(output_dir, 'LinDdef.pdf') if case == -1 else
                os.path.join(output_dir, f'NLD_{case}.pdf'))

    plt.savefig(filename, dpi=300, bbox_inches='tight')
    
    plt.close()

# Function: Clear Screen
def clear_screen():
    if platform.system() == "Windows":
        os.system('cls')
    else:
        os.system('clear')    
    
if __name__ == "__main__":
    clear_screen()
    main()



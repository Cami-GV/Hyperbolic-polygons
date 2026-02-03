'''
Welcome! Next you will find a Sage code for two main tasks:
-To define a GAP group with generators of prescribed order
-To draw a fundamental polygon for a surface endowed with a group action

Warnings:
-All the functions are defined on the unit disc. The metric on the disc is ds = 2|dz|/(1-|z|^2).
-The constructions are adapted for the case of groups acting with a signature of the form (0;m,n,2).

References:
-Anderson, J. Hyperbolic Geometry. Springer Undergraduate Mathematics Series, 2005
-Beardon, A. The Geometry of Discrete Groups. Graduate Texts in Mathematics, 91. Springer-Verlag, New York, 1983.
'''

from io import StringIO

z,w,vertex,c = var('z,w,vertex,c', domain='complex')
t,r,R,e = var('t,r,R,e', domain='real')
C = circle((0,0), 1, rgbcolor=(0,0,0))

Id = Matrix(CC,[[1,0],[0,1]])


#---Hyperbolic geometry---#

def apply_mobius(M: matrix, z:complex) -> complex:
    '''
    Input:
    -M: 2x2 complex matrix, representing a Möbius transformation
    -z: Complex number in the unit disc

    Output:
    -Mz: Image of z under M

    Example: 
        if M=[[a,b],[c,d]], apply_mobius(M,z)=(az+b)/cz+d
    '''
    Mz = (M[0][0]*z+M[0][1])/(M[1][0]*z+M[1][1])
    return Mz

def invert_mobius(M:matrix) -> matrix:
    '''
    Return a matrix N, representing the inverse Möbius transformation of M
    '''
    N = [[M[1][1],-M[0][1]],[-M[1][0],M[0][0]]]
    return N

def apply_translation(c:complex , z:complex) -> complex:
    '''
    A hyperbolic translation if a Möbius transformation of the form
    t -> (t-c)/(1-t*c.conjugate())

    Input:
    -c: Complex number. Represents the preimage of 0 
    -z: Complex number

    Output:
    -Mz: Image of z under a hyperbolic translation sending c to 0
    '''
    M = [[1,-c],[-c.conjugate(),1]] 
    Mz = apply_mobius(M,z)
    return Mz

def hyperbolic_length(z: complex, w:complex) -> float : 
    '''
    Return the hyperbolic distance between z and w

    Procedure:
    1. Map z to T(z)=0 and w to T(w) by a hyperbolic translation (isometry of the unit disc)
    2. Compute the distance between 0 and T(w)
    '''
    a = apply_translation(z,w)
    l = log((1+a.abs())/(1-a.abs())) 
    return l

def hyperbolic_angle(z: complex, w:complex, vertex:complex) -> float:
    '''
    Return the angle between the vectors z - vertex and w - vertex

    Procedure:
    1. Map vertex to T(vertex)=0 with a hyperbolic translation (isometry of the unit disc)
    2. Comput the angle between T(z) and T(w)
    '''
    Tz = apply_translation(vertex,z)
    Tw = apply_translation(vertex,w)
    angle = arccos((Tz.real()*Tw.real()+Tz.imag()*Tw.imag())/(Tz.abs()*Tw.abs()))
    return angle

def hyperbolic_geodesic(z:complex, w:complex, color=[]) -> plot:
    '''
    Input:
    -z: Complex number. It is the starting point of the geodesic
    -w: Complex number. It is the ending point of the geodesic
    -color: 

    Output:
    -P: Graphic object representing the hyperbolic geodesic from z to w

    Method:
        The geodesic starting at z and ending at w is the curve t -> (t*T(w)+z)/(1+t*z.conjugate()*T(w)),
        T is Möbius transformation sending z to 0
    '''
    T = apply_translation(z,w)
    P = parametric_plot([apply_translation(-z,t*T).real(), apply_translation(-z,t*T).imag()], (t,0,1), rgbcolor=color, axes=false)
    return P

def draw_Polygon(points:list, colors=[]) -> plot:
    '''
    Input:
    -points: List of complex numbers, representing the vertices of a polygon
    -colors: List of strings for color codes.
        By default, all edges are blue

    Output:
    -poly: Graphic object, representing a hyperbolic polygon with given vertices
    '''

    if colors==[]:
        colors = [[0,0,1] for k in range(len(points))]
    edges = [hyperbolic_geodesic(points[k],points[(k+1)%len(points)],colors[k]) for k in range(len(points))]
    poly = sum(edges)
    return poly

def rotate_Polygon(vertices: list, rotation: matrix) -> list:
    '''
    Input: 
    -vertices: List of complex numbers, representing the vertices of a polygon
    -rotation: Matrix, representing a Möbius transformation

    Output:
    -rotated_vertices: List of complex numbers. 
        Contains the images of vertices under the transformation rotation
    '''
    rotated_vertices = [apply_mobius(rotation,vertices[k]) for k in range(len(vertices))]
    return rotated_vertices

#---Fundamental polygons for signatures of type (0;m,n,2)---#
def vertices_Polygon_mn2(Sig: list) -> list:
    '''
    Input: 
    -Sig: -Sig: List of the form [m,n,2], which comes from a signature of type (0;m,n,2)

    Output: 
    -vertices: Four-complex numbers list, representing the vertices of a fundamental quadrilateral
    for the action with given signature 

    Description:
        The formulas are derived from the hyperbolic trigonometric laws.
        The desired fundamental quatrilateral has angles 2pi/m, pi/2, 2pi/n, and pi/2. 
    '''
    r = sqrt((cos(pi*Sig[1]^(-1))-sin(pi*Sig[0]^(-1)))/(cos(pi*Sig[1]^(-1))+sin(pi*Sig[0]^(-1)))).n()
    R = sqrt(cos(pi*(Sig[0]^(-1)+Sig[1]^(-1)))/cos(pi*(Sig[1]^(-1)-Sig[0]^(-1)))).n()
    z_0 = CDF(R*cos(pi*Sig[0]^(-1)),R*sin(pi*Sig[0]^(-1)))
    z_1 = CDF(r*cos(2*pi*Sig[0]^(-1)),r*sin(2*pi*Sig[0]^(-1)))
    vertices = [0,CDF(r,0),z_0,z_1]
    return vertices

def draw_Polygon_mn2(Sig: list,color=[]) -> list:
    '''
    Input: 
    -Sig: List of the form [m,n,2], which comes from a signature of type (0;m,n,2)
    -color: Four-string list for color codes.
        By default the colors are orchid-gray-gray-orchid

    Output: 
    -P: Graphic object. Rrepresents a fundamental quadrilateral for the action 
    on the unit disc with given signature.
    '''
    if len(color)==0:
        color = ['#DA70D6','#696969','#696969','#DA70D6']
    vertices = vertices_Polygon_mn2(Sig)
    P = draw_Polygon(vertices, color)
    return P

def rotations_mn2(Sig: list) -> list:
    '''
    Input: 
    -Sig: List of the form [m,n,2], which comes from a signature of type (0;m,n,2)

    Output: 
    - Two-matrices list, representing the generators of the uniformizing Fuchsian group
        for the quotient sphere, obtained under the action with given signature
    
    Description:
        The matrices are rotations of order m and n, centered at 0 and z_0 respectively 
        (see vertices_Polygon_mn2)
    '''
    vertices = vertices_Polygon_mn2(Sig) #relevant vertices: 0 and 2
    angles = [CDF(cos(2*pi*(Sig[k]^(-1))),sin(2*pi*(Sig[k]^(-1)))) for k in range(len(Sig))]
    A = Matrix(CC,[
        [angles[0],0],
        [0,1]
        ])
    B = Matrix(CC,[
        [angles[1]-(vertices[2].abs()^2), vertices[2]*(1-angles[1])],
        [-vertices[2].conjugate()*(1-angles[1]), 1-(angles[1]*(vertices[2].abs()^2))]
        ])
    return [A,B]

def fundamentalpolygon_mn2(id,Sig,colors=[]) -> plot:
    '''
    Input: 
    -id: GAP catalogue ID of a finite group
    -Sig: List of the form [m,n,2], which comes
        from a signature of type (0;m,n,2)
    -colors: Four-strings list for color codes.
        By default the colors are orchid-gray-gray-orchid
    
    Output: Graphic object, representing a fundamental polygon for an action
        of the group with given ID and signature in the unit disc.
    ''' 
    if len(colors) == 0:
        colors = ['#DA70D6','#696969','#696969','#DA70D6']
    K = with_signature_generators(id,Sig)
    matrices = as_matrices_mn2(K,Sig)
    vertices = vertices_Polygon_mn2(Sig)
    rotated_vertices = [rotate_Polygon(vertices,matrices[k]) for k in range(len(matrices))]
    rotated_pieces = [draw_Polygon(rotated_vertices[k], colors) for k in range(len(matrices))]
    poly = sum(rotated_pieces)
    return poly

#---Fundamental polygons signatures of type (0;m,2,2,2)---#

def vertices_m222(m: int, theta) -> list:
    '''
    Input:
    -m: Integer greater than 2
    -theta: an angle between 0 and pi, 
        representing an element in the real-parametrized family of the action. 

    Output:
    -vertices: Six-complex numbers list, representing the vertices of a fundamental hexagon
    for the action with signature (0;m,2,2,2). Two of the vertices are degenerate.

    Description:
        The formulas are derived from the hyperbolic trigonometric laws.
        The angles of the obtained polygon sum 2pi/m + 3pi.
    '''
    tau = pi/m
    r = sqrt((cos(tau)*tan(theta/2)-sin(tau)+1)/(cos(tau)*tan(theta/2)+sin(tau)+1)).n()
    R = sqrt((cos(tau)*cot(theta/2)-sin(tau)+1)/(cos(tau)*cot(theta/2)+sin(tau)+1)).n()
    cosht = sqrt((cos(tau)/sin(theta))+1)*sin(tau)^(-1)
    s = sqrt((cosht-1)/((cosht+1))).n()
    delta = acot((tan(theta/2)+cos(tau))/sin(tau))
    vertices = [0, r, 
                CDF(s*cos(delta),s*sin(delta)),
                CDF(R*cos(tau),R*sin(tau)),
                CDF(s*cos(2*tau-delta),s*sin(2*tau - delta)),
                CDF(r*cos(2*tau),r*sin(2*tau))]
    return vertices

def hexagon_m222(m: int, theta, colors=[]) -> plot:
    '''
    Input:
    -m: Integer, greater than 2
    -theta: an angle between 0 and pi, 
        representing an element in the real-parametrized family of the action
     -colors: Six-strings for color codes.
        By default the colors are orchid-silver-black-aqua-navy-orchid

    Output:
    -quad: Graphic object representing a fundamental polygon for an action with signature (0;m,2,2,2).
    '''
    if len(colors)==0:
        colors = ['orchid', 'silver', 'black', 'aqua', 'navy','orchid']
    verts = vertices_m222(m, theta)
    hex = draw_Polygon(verts, colors)
    return hex

def rotations_m222(m: int, theta) -> list:
    '''
    Input: 
    -m: Integer, grater than 2
    -theta: an angle between o and pi,
        representing and element in the real-parametrized family of the action.

    Output: Four-matrices list, representing the generators of a Fuchsian group 
        acting with signature (0;m,2,2,2)
    '''
    verts = vertices_m222(m,theta)
    X = Matrix(CC,[[CDF(cos(2*pi/m),sin(2*pi/m)),0],[0,1]]) #order m
    Y = Matrix(CC,[[-(1+verts[2].abs()^2), 2*verts[2]],[-2*verts[2].conjugate(), 1+verts[2].abs()^2]]) #order 2
    Z = Matrix(CC,[[-(1+verts[4].abs()^2), 2*verts[4]],[-2*verts[4].conjugate(), 1+verts[4].abs()^2]]) #order 2
    W = X*Y*Z
    return [X,Y,Z,W]

def fundamental_polygon_m222(id, Sig, theta, colors=[]) -> plot:
    '''
    Input: 
    -id: GAP catalogue ID for a finite group
    -Sig: List of the form [m,2,2,2], which comes from a signature
        of type (0;m,2,2,2)
    -theta: Angle between 0 and pi, representing a member for the real-parametrized
        family with the given action
    -colors: Six-strings for color codes.
        By default the colors are orchid-silver-black-aqua-navy-orchid

    Output: Graphic object, representing a fundamental polygon for an action
        of the group with given ID and signature in the unit disc.
    '''
    if len(colors)==0:
        colors = ['orchid', 'silver', 'black', 'aqua', 'navy','orchid']
    K = with_signature_generators(id,Sig)
    matrices = as_matrices_m222(K,Sig,theta)
    verts = vertices_m222(Sig[0], theta)
    rotated_vertices = [rotate_Polygon(verts,matrices[k]) for k in range(len(matrices))]
    rotated_pieces = [draw_Polygon(rotated_vertices[k],colors) for k in range(len(rotated_vertices))]
    poly = sum(rotated_pieces)
    return poly

#---Groups manipulations---#

def fp_representation(id):
    '''
    Input: 
    -id: GAP catalogue ID of a finite group

    Output: 
    -G2: Finitely presented GAP group, which is isomorphic  
        to the group with the given ID
    '''
    G = libgap.SmallGroup(id[0],id[1])
    G1 = libgap.IsomorphismFpGroup(G).Image()
    G2 = libgap.IsomorphismSimplifiedFpGroup(G1).Image()
    return G2

def search_generators(G,Sig:list) -> list:
    '''
    Input: 
    -G: Finitely presented GAP group
    -Sig: List of the form [m,n,2], which comes
        from a signature of type (0;m,n,2)

    Output: 
    -search: List of lists. 
        The i-th list contains all elements of G of order Sig[i].
    '''
    group_elements = libgap.Elements(G)
    search = []
    for m in range(len(Sig)):
        could_be = []
        for k in range(len(group_elements)):
            if libgap.Order(group_elements[k])==Sig[m]:
                could_be.append(group_elements[k])
        search.append(could_be)
    return search

generators = []
def potencial_generators(G,search:list, acumulados:list=[]) -> list:
    '''
    Input: 
    -G: Finitely presented GAP group
    -search: List obtained from search_generators.
        Contains elements of G with orders prescribed by a signature
    -acumulados: variable list. By defaul is the empty list 

    Output: 
    -generators: List of all possible sets of generators of G 
    with the prescribed orders
    '''
    if len(search)==0:
        if prod(acumulados) == libgap.One(G):
            S = libgap.Subgroup(G, acumulados)
            if libgap.Index(G,S)==1:
                generators.append(acumulados)
    else:
        for g in search[0]:
            potencial_generators(G, search[1:], acumulados+[g])
    return generators

def with_signature_generators(id, Sig:list):
    '''
    Input: 
    -id: GAP catalogue ID of a finite group
    -Sig: List of the form [m,n,2], which comes from
        a signature of type (0;m,n,2)

    Output: 
    -K: Finitely presented GAP group isomorphic to the given ID,
        and with prescribed generators
    '''
    G = fp_representation(id)
    generators.clear()
    search = search_generators(G,Sig)
    gens = potencial_generators(G,search)
    K = libgap.IsomorphismFpGroupByGenerators(G,gens[0]).Image()
    return K

def as_matrices_mn2(K,Sig) -> list:
    '''
    Input:
    -K: Finitely presented GAP group, obtained using With_signature_generators 
    -Sig: List of the form [m,n,2], which comes from
        a signature of type (0;m,n,2)

    Output: 
    -elements_as_mat: List containing matrices.
        These matrices represent the elements of K as Möbius transformations 
    '''
    s = StringIO()
    elements = libgap.Elements(K)
    print(elements,file=s)
    rotations = rotations_mn2(Sig)
    X = rotations[0]
    Y = rotations[1]
    Z = rotations[0]*rotations[1]
    elements_str = [f.strip() for f in s.getvalue()[1:-2].replace('F1', 'X').replace('F2', 'Y').replace('F3','Z').replace('^', '**').split(',')]
    elements_as_mat = [Id]+[sage_eval(expr, locals = {'X':X, 'Y':Y, 'Z':Z}) for expr in elements_str[1:]]
    return elements_as_mat

def as_matrices_m222(K,Sig,theta) -> list:
    '''
    Input:
    -K: Finitely presented GAP group, obtained using With_signature_generators 
    -Sig: List of the form [m,2,2,2], which comes from
        a signature of type (0;m,2,2,2)

    Output: 
    -elements_as_mat: List containing matrices.
        These matrices represent the elements of K as Möbius transformations 
    '''
    s = StringIO()
    elements = libgap.Elements(K)
    print(elements,file=s)
    rotations = rotations_m222(Sig[0],theta)
    X = rotations[0]
    Y = rotations[1]
    Z = rotations[2]
    W = rotations[3]
    elements_str = [f.strip() for f in s.getvalue()[1:-2].replace('F1', 'X').replace('F2', 'Y').replace('F3','Z').replace('F4','W').replace('^', '**').split(',')]
    elements_as_mat = [Id]+[sage_eval(expr, locals = {'X':X, 'Y':Y, 'Z':Z, 'W':W}) for expr in elements_str[1:]]
    return elements_as_mat

def as_word_in_G(K,elem):
    '''
    Input:
    -K: Finitely presented GAP group
    -elem: Expression in the generators of K.

    Ouput:
    - An element of K, equivalent to elem under the relations of the group
    '''
    element = libgap.Elements(K)
    factor = libgap.Factorization(K,elem) #returns an abstract expression in the generators of K, not a GAP element
    for k in range(len(element)):
        if libgap.Factorization(K,element[k])==factor:
            return element[k]

#Specifically for Bolza surface
S = Matrix(ZZ, [[0,0,0,-1],[1,0,0,0],[0,1,0,0],[0,0,1,0]]) #x action
T = Matrix(ZZ,[[0,-1,0,0],[1,-1,0,0],[1,0,-1,1],[0,1,-1,0]]) #y action

def action_of_G_mn2(K,actions) -> list:
    '''
    Input:
    -K: Finitely presented GAP group
    -actions: Two-matrices list.
        Contains the representing matrices of the action of the generators of K
    
    Output:
    -actions_as_mat: List of matrices.
        Contains the representing matrices for the actions of the elements of K
    
    '''
    w= StringIO()
    elements = libgap.Elements(K)
    print(elements,file=w)
    S = actions[0]
    T = actions[1]
    U = S*T
    actions_str = [f.strip() for f in w.getvalue()[1:-2].replace('F1', 'S').replace('F2', 'T').replace('F3','U').replace('^', '**').split(',')]
    actions_as_mat = [Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])]+[sage_eval(expr, locals = {'S':S, 'T':T, 'U':U}) for expr in actions_str[1:]]
    return actions_as_mat

def stabilizers_mn2(K,vec,actions) -> list:
    '''
    Input:
    -K: Finitely presented GAP group
    -vec: Four-integer vector
    -actions: List of matrices.
        Constains the matrices representing the action of the generators of K 
        on a given homology basis

    Output:
    -stab_elements: List of matrices.
        Contain the elements of K fixing vec under its corresponding action
    '''
    elem = libgap.Elements(K)
    stabs = []
    Actions = action_of_G_mn2(K,actions)
    for k in range(len(Actions)):
        if Actions[k]*vec==vec or Actions[k]*vec==-vec:
            stabs.append([Actions[k],k])
    print(len(stabs))
    stab_elements = []
    for m in range(len(stabs)):
        M = stabs[m]
        stab_elements.append(elem[M[1]])
    return stab_elements

def action_of_G_m222(K,actions) -> list:
    '''
    Input:
    -K: Finitely presented GAP group
    -actions: Two-matrices list.
        Contains the representing matrices of the action of the generators of K
    
    Output:
    -actions_as_mat: List of matrices.
        Contains the representing matrices for the actions of the elements of K
    
    '''
    w= StringIO()
    elements = libgap.Elements(K)
    print(elements,file=w)
    S = actions[0]
    T = actions[1]
    U = actions[2]
    V = S*T*U
    actions_str = [f.strip() for f in w.getvalue()[1:-2].replace('F1', 'S').replace('F2', 'T').replace('F3','U').replace('F4','V').replace('^', '**').split(',')]
    actions_as_mat = [Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])]+[sage_eval(expr, locals = {'S':S, 'T':T, 'U':U, 'V':V}) for expr in actions_str[1:]]
    return actions_as_mat

def stabilizers_m222(K,vec,actions) -> list:
    '''
    Input:
    -K: Finitely presented GAP group
    -vec: Four-integer vector
    -actions: List of matrices.
        Constains the matrices representing the action of the generators of K 
        on a given homology basis

    Output:
    -stab_elements: List of matrices.
        Contain the elements of K fixing vec under its corresponding action
    '''
    elem = libgap.Elements(K)
    stabs = []
    Actions = action_of_G_m222(K,actions)
    for k in range(len(Actions)):
        if Actions[k]*vec==vec or Actions[k]*vec==-vec:
            stabs.append([Actions[k],k])
    print(len(stabs))
    stab_elements = []
    for m in range(len(stabs)):
        M = stabs[m]
        stab_elements.append(elem[M[1]])
    return stab_elements


#---Classes---#

class SignatureGroup():
    '''
    Class of GAP groups with prescribed generators by signature

    Methods:
    -ID: GAP catalogue ID for a finite group
    -signature: List of the form [m,n,p], which comes from a signature of type (d; m,n,p)
    -signature_generated: Finitely presented GAP group with generators prescribed by signature
    -elements: List of GAP elements of the group
    -generators: List of generators of the GAP group
    -relations: List of relations of the finitely presented GAP group
    
    -reduction: Perform reductions to obtain an equivalent expression of a given word,
        expressed in the generators of the group
    '''
    def __init__(self, id, Sig):
        self.ID = id
        self.order = id[0]
        self.signature = Sig
        self.signature_generated = with_signature_generators(self.ID, self.signature)
        self.elements = libgap.Elements(self.signature_generated)
        self.generators = libgap.GeneratorsOfGroup(self.signature_generated)
        self.relations = libgap.RelatorsOfFpGroup(self.signature_generated)

    def reduction(self, word):
        '''
        Returns an element in the group equivalent to the given word.

        Input: a GAP word, expressed in the generators of the group

        Example:
        sage: G = SignatureGroup([48,29],[8,3,2])
        sage: gens = G.generators
        sage: x = gens[0]
        sage: y = gens[1]
        sage: G.reduction(y^2)
        F1*F2*F1
        '''
        factor = libgap.Factorization(self.signature_generated, word)
        for k in range(len(self.elements)):
            if libgap.Factorization(self.signature_generated, self.elements[k]) == factor:
                return self.elements[k]

class FundamentalPolygon_mn2():
    '''
    Class of fundamental polygons for a given group and signature

    Methods:
    -ID: GAP catalogue ID for a finite group
    -signature: List of the form [m,n,2], which comes form a signature of type (0;m,n,2)
    
    -group: Finitely presented GAP group with generators prescribed by signature
    
    -quadrilateral_vertices: Four-complex numbers list. 
        Contains the vertices of a fundamental quadrilateral for the quotient sphere
    -color: Four.strings list for color codes. 
    -rotations_quadrilateral: Two-matrices list, representing, as Möbius transformations, 
        the generators of the uniformizing Fuchsian group for the quotient sphere
    -G_as_matrices: List of matrices.
        Contains a matrix representatives for the action in homology of each element 
        of the group acting on the surface
    
    -genus: Return the genus on which the given group acts with the given signature
    -draw_quadrilateral: Return a graphic object, representing the fundamental quadrilateral for the quotient sphere
    -draw: Return a graphic object, which represents a fundamental polygon for a surface endowed with the given group action
    -G_action: Return a list of 2gx2g matrices, where g is the genus of the surface.
        Contains the matrix representative fo the action in homology of each element of the group
    '''
    def __init__(self, id, Sig):
        self.ID = id
        self.signature = Sig 
        self.group = SignatureGroup(id, Sig)
        self.quadrilateral_vertices = vertices_Polygon_mn2(Sig)
        self.color = ["orchid", "#696969", "#696969", "orchid"]
        
        self.rotations_quadrilateral = rotations_mn2(self.signature)
        self.X = self.rotations_quadrilateral[0]
        self.Y = self.rotations_quadrilateral[1]

        self.G_as_matrices = as_matrices_mn2(self.group.signature_generated, self.signature) 

    def change_color(self, colors):
        self.color = colors

    def genus(self):
        #Computed using Riemann-Hurwitz formula
        '''
        Example:
        sage: Poly = FundamentalPolygon_mn2([48,29],[8,3,2])
        sage: Poly.genus()
        2
        '''
        order = self.group.order
        len_sig = len(self.signature)
        sig_sum = sum([self.signature[k]^-1 for k in range(len(self.signature))])
        g = 1 + order*((2^-1)*(len_sig-sig_sum)-1)
        return g

    def draw_quadrilateral(self):
        '''
        Example:
        sage: Poly = FundamentalPolygon_mn2([48,29],[8,3,2])
        sage: Poly.draw_quadrilateral()
        '''
        return draw_Polygon_mn2(self.signature)
    
    def draw(self, colors: list=[]) -> plot:
        '''
        Input:
        -colors: Four-string list of color codes. 
            By default the colors are orchid-gray-gray-orchid
        Output:
        -Poly: Graphic object, representing a fundamental polygon for an action
        of the group with given ID and signature in the unit disc.

        Example:
        sage: Poly = FundamentalPolygon_mn2([48,29],[8,3,2])
        sage: Poly.draw()
        ''' 
        if len(colors)==0:
            colors = ['#DA70D6','#696969','#696969','#DA70D6']
        rotated_vertices = [rotate_Polygon(self.quadrilateral_vertices, self.G_as_matrices[k]) for k in range(len(self.G_as_matrices))]
        rotated_pieces = [draw_Polygon(rotated_vertices[k], colors) for k in range(len(rotated_vertices))]
        Poly = sum(rotated_pieces)
        return Poly

    def G_action(self, actions: list):
        '''
        Input:
        -actions: List of two 2gx2g matrices, representing the action of the generators of the group 
            in a given homology basis.
        
        Output:
        -actions_as_mat: List of matrices, containing the representations of the action in homology 
            of a each element of the group.
        '''
        w = StringIO()
        print(self.group.elements, file=w)
        S = actions[0]
        T = actions[1]
        U = S*T
        actions_str = [f.strip() for f in w.getvalue()[1:-2].replace('F1', 'S').replace('F2', 'T').replace('F3','U').replace('^', '**').split(',')]
        actions_as_mat = [Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])]+[sage_eval(expr, locals = {'S':S, 'T':T, 'U':U}) for expr in actions_str[1:]]
        return actions_as_mat

class FundamentalPolygon_m222():
    '''
    Class of fundamental polygons for a given group and signature

    Methods:
    -ID: GAP catalogue ID for a finite group
    -signature: List of the form [m,n,2], which comes form a signature of type (0;m,n,2)
    
    -group: Finitely presented GAP group with generators prescribed by signature
    
    -quadrilateral_vertices: Four-complex numbers list. 
        Contains the vertices of a fundamental quadrilateral for the quotient sphere
    -color: Four.strings list for color codes. 
    -rotations_quadrilateral: Two-matrices list, representing, as Möbius transformations, 
        the generators of the uniformizing Fuchsian group for the quotient sphere
    -G_as_matrices: List of matrices.
        Contains a matrix representatives for the action in homology of each element 
        of the group acting on the surface
    
    -genus: Return the genus on which the given group acts with the given signature
    -draw_quadrilateral: Return a graphic object, representing the fundamental quadrilateral for the quotient sphere
    -draw: Return a graphic object, which represents a fundamental polygon for a surface endowed with the given group action
    -G_action: Return a list of 2gx2g matrices, where g is the genus of the surface.
        Contains the matrix representative fo the action in homology of each element of the group
    '''
    def __init__(self, id, Sig, theta):
        self.ID = id
        self.m = Sig[0] 
        self.signature = Sig
        self.parameter = theta
        self.group = SignatureGroup(id, Sig)
        self.hexagon_vertices = vertices_m222(self.m, self.parameter)
        self.color = ['orchid','silver', 'black','aqua','navy','orchid']
        
        self.rotations_hexagon = rotations_m222(self.m,self.parameter)
        self.X = self.rotations_hexagon[0]
        self.Y = self.rotations_hexagon[1]
        self.Z = self.rotations_hexagon[2]

        self.G_as_matrices = as_matrices_m222(self.group.signature_generated, self.signature, self.parameter) 

    def change_color(self, colors):
        self.color = colors

    def genus(self):
        #Computed using Riemann-Hurwitz formula
        '''
        Example:
        sage: Poly = FundamentalPolygon_mn2([12,4],[3,2,2,2])
        sage: Poly.genus()
        2
        '''
        order = self.group.order
        len_sig = len(self.signature)
        sig_sum = sum([self.signature[k]^-1 for k in range(len(self.signature))])
        g = 1 + order*((2^-1)*(len_sig-sig_sum)-1)
        return g

    def draw_hexagon(self):
        '''
        Example:
        sage: Poly = FundamentalPolygon_mn2([48,29],[8,3,2])
        sage: Poly.draw_quadrilateral()
        '''
        return hexagon_m222(self.m, self.parameter)
    
    def draw(self, colors: list=[]) -> plot:
        '''
        Input:
        -colors: Four-string list of color codes. 
            By default the colors are orchid-gray-gray-orchid
        Output:
        -Poly: Graphic object, representing a fundamental polygon for an action
        of the group with given ID and signature in the unit disc.

        Example:
        sage: Poly = FundamentalPolygon_mn2([48,29],[8,3,2])
        sage: Poly.draw()
        ''' 
        if len(colors)==0:
            colors = ['orchid','silver','black','aqua','navy','orchid']
        rotated_vertices = [rotate_Polygon(self.hexagon_vertices, self.G_as_matrices[k]) for k in range(len(self.G_as_matrices))]
        rotated_pieces = [draw_Polygon(rotated_vertices[k], colors) for k in range(len(rotated_vertices))]
        Poly = sum(rotated_pieces)
        return Poly

    def G_action(self, actions: list):
        '''
        Input:
        -actions: List of two 2gx2g matrices, representing the action of the generators of the group 
            in a given homology basis.
        
        Output:
        -actions_as_mat: List of matrices, containing the representations of the action in homology 
            of a each element of the group.
        '''
        w = StringIO()
        print(self.group.elements, file=w)
        S = actions[0]
        T = actions[1]
        U = actions[2]
        V = S*T*U
        actions_str = [f.strip() for f in w.getvalue()[1:-2].replace('F1', 'S').replace('F2', 'T').replace('F3','U').replace('F4','V').replace('^', '**').split(',')]
        actions_as_mat = [Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])]+[sage_eval(expr, locals = {'S':S, 'T':T, 'U':U, 'V':V}) for expr in actions_str[1:]]
        return actions_as_mat



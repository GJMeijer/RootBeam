# Geometry functions - version 01/09/2017

#IMPORT PACKAGES
import numpy as np
import module_Geometry as mGeom
import module_SoilRoot as mSoil


#LOCAL TO GLOBAL SEGMENT POSITIONS (dimensional)
def func_localtoglobalsegmentpositions(x,y,X1,Y1,Theta):
    """
    FUNCTION to calculate (X,Y) (GLOBAL positions) from local dimensional positions (x,y)
    INPUT
    - <x>:     Local x-coordinate (aligns with undisplaced root axis), dimensional, x=0 at Node1
    - <y>:     Local y-coordinate (perpendicular to  undisplaced root axis), dimensional
    - <X1>:    Global X-coordinate of segment Node1
    - <Y1>:    Global Y-coordinate of segment Node1
    - <Theta>: Orientation of undisplaced root axis in Global X-Y system
    OUTPUT
    - <X, Y>:  arrays with Global positions
    """
    X = X1 + x * np.cos(Theta) - y * np.sin(Theta)
    Y = Y1 + x * np.sin(Theta) + y * np.cos(Theta)
    return(X, Y)


#LOCAL TO GLOBAL DISPLACEMENTS (dimensional)
def func_localdeformation2globaldeformations(u,w,Theta):
    """
    FUNCTION to calculate (u,w) (LOCAL displacements) from GLOBAL displacements (U,W)
    INPUT
    - <u>:     Local displacement in x-direction (dimensinal)
    - <w>:     Local displacement in y-direction (dimensional)
    - <Theta>: Orientation of undisplaced root axis in Global X-Y system
    OUTPUT
    - <U, W>:  arrays with dimensional displacements in Global coordinate system
    """
    U = u * np.cos(Theta) - w * np.sin(Theta)
    W = u * np.sin(Theta) + w * np.cos(Theta)
    return(U, W)   


#GLOBAL TO LOCAL DISPLACEMENTS (dimensional)
def func_globaldeformation2localdeformations(U,W,Theta):
    """
    FUNCTION to calculate (u,w) (LOCAL displacements) from GLOBAL displacements (U,W)
    INPUT
    - <U>:     Global displacement in X-direction
    - <W>:     Global displacement in Y-direction
    - <Theta>: Orientation of undisplaced root axis in Global X-Y system
    OUTPUT
    - <u, w>:  arrays with dimensional displacements in local coordinate system
    """
    u = U *  np.cos(Theta) + W * np.sin(Theta)
    w = U * -np.sin(Theta) + W * np.cos(Theta)
    return(u, w)   


#CALCULATE GLOBAL POSITIONS BASED ON LOCAL, DIMENSIONLESS DISPLACEMENTS FOR ALL SEGMENTS
def func_localdimensionless2globalpositions(v, dsegm):
    """
    FUNCTION to calculate global positions based on dimensionless local positions
    INPUT
    - <v>:     Current displacement state of root (local coordinate system, dimensional parameters)
    - <dsegm>: Segment data
    OUTPUT
    - <out>:   dictionary with global positions (Xr and Yr) for every segment
    """
    #output dictionary
    out = dict()
    #loop through segments
    for sID,i in v.items():
        #per segment, create dictionary
        out[sID] = dict()
        #local positions - dimensionless
        x = i['u'] + i['s']
        y = i['w']
        #global positions
        Xr, Yr = func_localtoglobalsegmentpositions(x, y, dsegm[sID]['X1'], dsegm[sID]['Y1'], dsegm[sID]['Theta'])
        out[sID]['Xr'] = Xr
        out[sID]['Yr'] = Yr 
       
    #return
    return(out)


#INITIAL ROOT POSITIONS
def func_initialglobalrootposition(dsegm, p):
    """
    FUNCTION to calculate initial root positions
    INPUT
    - <dsegm>:  segment data
    - <p>:      Model input parameters
    OUTPUT
    - <out>:    dictionary with intial global root position
    """
    #output dictionary
    out = dict()
    #loop through segments
    for sID,i in dsegm.items():
        out[sID] = dict()
        out[sID]['Xr0'], out[sID]['Yr0'] = func_localtoglobalsegmentpositions(np.linspace(0,i['L'],p['n_node']), np.linspace(0,0,p['n_node']), i['X1'], i['Y1'], i['Theta'])
    #return
    return(out)

#DEFORMATIONS AT NODE
def func_nodaldeformations(uext, dnode, p):
    """
    FUNCTION to get soil deformations at original positions of node
    INPUT
    - <uext>:   shear displacement
	- <dnode>:  dictionary with nodal data
    - <p>:      Model input parameters
    OUTPUT
    - <out>:    dictionary with soil deformations at original positions of node
    """
    out = dict()
    #create dictionary entries
    for nID,i in dnode.items():
        out[nID] = dict()
    #get soil deformations at original nodal positions
    for nID,i in out.items():
        ush,wsh = mSoil.func_soildisplacement_shear(dnode[nID]['X'], dnode[nID]['Y'], uext, p)
        i['Us'],i['Ws'] = mGeom.func_localdeformation2globaldeformations(ush, wsh, p['shearplane_angle'])
    return(out)   


#DOUBLE ASYMPTOTE FUNCTIONS - SMOOTHING
def zeta1(x, b, xd=0.0):
    """
    FUNCTION for double asymptote
    - x -> -inf,  y -> 0
    - x ->  inf   y -> x - xd
    - deflection points at x=xd
    INPUT
    - <x>:     x-coordinate (numpy array)
    - <b>:     parameter indicating the rapidness of change of curve
    OPTIONAL INPUT
    - <xd>:    x-coordinate of mid-point where curve changes direction
    OUTPUT
    - <y>:     y-values
    """
    return(0.5 * ((x-xd) + (1.0+(b*(x-xd))**2)**0.5 / b))

def zeta2(x, b, xd=0.0):
    """
    FUNCTION for double asymptote:
    zeta2 = d/dx(zeta1)
    - x -> -inf,  y -> 0
    - x ->  inf   y -> 1
    - deflection points at x=xd
    INPUT
    - <x>:     x-coordinate (numpy array)
    - <b>:     parameter indicating the rapidness of change of curve
    OPTIONAL INPUT
    - <xd>:    x-coordinate of mid-point where curve changes direction
    OUTPUT
    - <y>:     y-values
    """
    return(0.5 + 0.5*b*(x-xd) * (1.0+(b*(x-xd))**2)**-0.5)

def zeta3(x, b, xd=0.0):
    """
    FUNCTION for double asymptote:
    zeta3 = d/dx(zeta2) = d^2/dx^2(zeta1)
    - x -> -inf,  y -> 0
    - x ->  inf   y -> 0
    - deflection points at x=xd
    INPUT
    - <x>:     x-coordinate (numpy array)
    - <b>:     parameter indicating the rapidness of change of curve
    OPTIONAL INPUT
    - <xd>:    x-coordinate of mid-point where curve changes direction
    OUTPUT
    - <y>:     y-values
    """
    return(0.5*b * (1.0+(b*(x-xd))**2)**-1.5)
# Soil and root functions - version 19/09/2017
# takes effect of increased lateral stress on axial effects into account
# takes effect of root rotation into account

#IMPORT PACKAGES
import numpy as np
import module_Geometry as mGeom


#FUNCTION TO GET CURRENT ROOT SECANT STIFFNESS
def func_secantstiffness(eps, E1, E2, E3):
    """
    FUNCTION to calculate secant stiffness after strain eps
    INPUT
    - <eps>:         root strain (for bending: strain in ultimate fibre)
    - <E1, E2, E3>:  fitting parameters
    OUTPUT
    - <Esec>:        Secant stiffness
    
    NOTES:
    - stress is fitted in the form:
          sig = E1*eps / sqrt(1+E2*eps^2) + E3*eps
      therefore:
          Esec = sig / eps = E1 / sqrt(1+E2*eps^2) + E3
    - This formulation for stress has the benefit that 
      1: the secant stiffness is finite for eps=0
      2: correct stiffnesses are calculated for negative strains
    """
    
    #calculate secant stiffness
    Esec = E3 + E1 / np.sqrt(1.0 + E2 * eps**2)
    #return
    return(Esec)


#FUNCTION TO CALCULATE GLOBAL SOIL DISPLACEMENT INCREMENT
def func_soildisplacement_shear(X, Y, uext, p):
    """
    FUNCTION to calculate U (GLOBAL X-direction) and W (GLOBAL Y-direction) of soil displacement at any global position (X,Y)
    INPUT
    - <X>: global X-coordinate where displacement is to be calculated (non-normalised)
    - <Y>: global Y-coordinate where displacement is to be calculated (non-normalised)
    - <uext>: magnitude of soil displacement (non-normalised)
    - <p>: dictionary with input parameters
    OUTPUT
    - <U, W>: arrays with displacement components in global X and Y-direction
    """
    
    #get shear plane depth at location (X,Y)
    Xs = p['shearplane_depth'] + Y * np.cos(p['shearplane_angle'])
    #distance to the shear plane (perpendicular, below the shearplane = positive)
    dist = (X - Xs) * np.sin(p['shearplane_angle']) - Y * np.cos(p['shearplane_angle'])
    #get the total shear displacement at location (X,Y) in direction of shearplane, using a hyperbolic tangent curve
    ush = uext * 0.5 * (np.tanh(dist / p['shearplane_param']))
    #displacement perpendicular to shear plane
    wsh = 0.0
    #return
    return(ush, wsh)


#MOBILISATION DISTANCE
def func_mobilisationdistance(uext, X, Y, u, w, t, dsegmsingle, p):
    """
    FUNCTION to calculate axial and lateral mobilisation distances
    INPUT
    - <u>:            root axial displacement (local coordinate system x-y, dimensional)
    - <w>:            soil lateral displacement (local coordinate system x-y, dimensional)
    - <dsegmsingle>:  dictionary with segment data for one single segment
    - <p>:            input parameter dictionary
    OUTPUT
    - <umob>, <wmob>: axial and lateral mobilisation distances (dimensional)
    
    CONCEPT
    - Axial mobilisation  - defined as change in distance to the shear plane
    - Lateral moblisation - defined as change in direction perpendicular to non-displaced root segment axis
    """
    #get soil displacements (in direction of shear plane
    ushear,wshear = func_soildisplacement_shear(X, Y, uext, p)
    #soil displacement in coordinate system aligned with deformed root    
    us2,ws2 = mGeom.func_localdeformation2globaldeformations(ushear, wshear, p['shearplane_angle']-dsegmsingle['Theta']-t)
    #soil displacements
    wtemp = w / np.sin(p['shearplane_angle'])
    u2    = wtemp * np.cos(p['shearplane_angle']-dsegmsingle['Theta']-t) 
    w2    = wtemp * np.sin(p['shearplane_angle']-dsegmsingle['Theta']-t)
    u3    = u - w * np.cos(p['shearplane_angle']-dsegmsingle['Theta'])
    #total mobilisation
    umob  = us2 - u2 - u3
    wmob  = ws2 - w2
    #return
    return(umob, wmob)


#FUNCTION TO CALCULATE CURRENT MOBILISATION PARAMETER
def func_mobilisationparameters(umob, wmob, ba, bl):
    """
    FUNCTION to calculate axial and lateral mobilisation fraction
    INPUT
    - <u>:              axial displacement (local coordinate system x-y, dimensional)
    - <w>:              lateral displacement (local coordinate system x-y, dimensional)
    - <dsegmsingle>:    dictionary with segment data for one single segment
    - <p>:              input parameter dictionary
    OUTPUT
    - <zetaa>, <zetal>: axial and lateral mobilisation fraction
    """
    #get mobilisation fractions
    zetaa  = np.tanh(umob / ba)
    zetal  = np.tanh(wmob / bl)
    #return
    return(zetaa, zetal)


#ABOVE SHEARPLANE CORRECTION
def func_mobilisationshearplanecorrection(X, Y, p):
    """
    FUNCTION to calculate axial and lateral mobilisation fraction
    INPUT
    - <x>:              current axial position (local coordinate system x-y, dimensional)
    - <y>:              current lateral position (local coordinate system x-y, dimensional)
    - <dnode>:          nodal data dictionary
    - <dsegmsingle>:    dictionary with segment data for one single segment
    - <p>:              input parameter dictionary
    OUTPUT
    - <zetash>:         mobilisation factor for position with respect to shearplane

    CONCEPT
    - If root crosses the shearplane, mobilisation direction is reversed, both for lateral and axial resistance
    """
    #Shear plane position
    Xsh = p['shearplane_depth']
    Ysh = 0.0
    #distance to shearplane
    xd, yd = mGeom.func_globaldeformation2localdeformations(X-Xsh, Y-Ysh, p['shearplane_angle']-0.5*np.pi)
    #correction factor - axial and lateral - resistance reduced when near shear plane
    zetasha = np.abs(np.tanh(xd / p['bsh']))
    zetashl = np.abs(np.tanh(xd / p['bsh']))
    #return
    return(zetasha, zetashl)
    

#CURRENT GLOBAL SOIL POSITIONS
def func_currentglobalsoilposition(V0, uext, p):
    """
    FUNCTION to calculate initial root positions
    INPUT
    - <V0>:     dictionary with initial position (global, dimensional)
    - <uext>:   external soil displacement
    - <p>:      Model input parameters
    OUTPUT
    - <out>:    dictionary with intial global root position
    """
    #output dictionary
    out = dict()
    #loop through segments
    for sID,i in V0.items():
        #dictionary with results
        out[sID] = dict()
        #displacement in direction of shearplane
        ush, wsh = func_soildisplacement_shear(i['Xr0'],  i['Yr0'], uext, p)
        #in global coordinate system
        Us, Ws   = mGeom.func_localdeformation2globaldeformations(ush, wsh, p['shearplane_angle'])
        #global root positions
        out[sID]['Xs'] = i['Xr0'] + Us
        out[sID]['Ys'] = i['Yr0'] + Ws
    #return
    return(out)


#CALCULATE CURRENT EFFECTIVE STRESS
def func_effectivestress(X, p, bx=1.0):
    """
    FUNCTION to calculate vertical effective soil stress at depth X
    INPUT
    - <X>       Depth (dimensional)
    - <p>:      Model input parameters
    OPTIONAL INPUT
    - <bx>:     curve parameter near surface, to smooth out sudden changes in stress (zeta functions)
    OUTPUT
    - <sigmav>: Soil vertical effective stress
    """
    return(mGeom.zeta1(X,bx) * p['gamma'] + mGeom.zeta2(X,bx) * p['sigmav0'])
    

#CALCULATE CURRENT EFFECTIVE STRESS GRADIENT
def func_effectivestressgradient(X, p, bx=1.0):
    """
    FUNCTION to calculate gradient in vertical effective soil stress at depth X
    INPUT
    - <X>       Depth (dimensional)
    - <p>:      Model input parameters
    OPTIONAL INPUT
    - <bx>:     curve parameter near surface, to smooth out sudden changes in stress (zeta functions)
    OUTPUT
    - <sigmav>: Soil vertical effective stress
    """
    return(mGeom.zeta2(X,bx) * p['gamma'] + mGeom.zeta3(X,bx) * p['sigmav0'])    


#Function to calulcate the mobilisation parameter for soil resistance
def func_soilresistance(uext, s, u, w, t, dsegmsingle, p):
    """
    FUNCTION to calculate soil resistance increments
    INPUT
    - <Xr>:          current position - depth, in global coordinate system
    - <dsegmsingle>: segment data
    - <delta_u>:     Mobilisation distance axial resistance   (dimensional, local coordinate system)
    - <delta_w>:     Mobilisation distance lateral resistance (dimensional, local coordinate system)
    - <p>:           Model input parameters
    OUTPUT
    - <qa>: Axial resistance increment
    - <ql>: Lateral resistance increment
    """
    
    #local positions
    x = s + u
    y = w
    #global positions
    X,Y = mGeom.func_localtoglobalsegmentpositions(x,y,dsegmsingle['X1'],dsegmsingle['Y1'],dsegmsingle['Theta'])
    
    #current vertical stress
    sigmav  = func_effectivestress(X, p, bx=2.0)

    #LATERAL RESISTANCE - REESE AND VAN IMPE
    #Reese&VanImpe parameters
    phi   = p['phi'] / 180*np.pi
    alpha = phi / 2.0
    beta  = np.pi/4.0 + alpha
    K0    = 0.4
    Ka    = (np.tan(np.pi/4.0 - alpha))**2
    
    #Equivalent depth
    Xreq = mGeom.zeta1(X, 1.0, xd=0.0) #+ p['sigmav0'] / p['gamma']
    pst = sigmav * (K0*Xreq*np.tan(phi)*np.sin(beta) / (np.tan(beta-phi)*np.cos(alpha)) + 
                    np.tan(beta)/np.tan(beta-phi) * (dsegmsingle['d']+Xreq*np.tan(beta)*np.tan(alpha)) + 
                    K0*Xreq*np.tan(beta) * (np.tan(phi)*np.sin(beta)-np.tan(alpha)) - 
                    Ka*dsegmsingle['d'])
    psd = sigmav * dsegmsingle['d'] * (Ka*(np.tan(beta)**8 - 1.0) + K0*np.tan(phi)*np.tan(beta)**4)
    #pd, equal to min(pst, psd) - per unit pile length
    pd = np.minimum(pst, psd)

    #As - static multiplication factor - some smoothing is required to get a right curve, hence the use of zeta1-functions
    As = (X/dsegmsingle['d'] - 4.5)**2 * 0.0988 + 0.88
    As[X < 0.0] = 4.5**2*0.0988 + 0.88
    As[X/dsegmsingle['d'] > 4.5] = 0.88

    #ultimate resistance according to reese&vanImpe - per unit length and diameter
    pu = As * pd / dsegmsingle['d']
    
    #MOBILISATION
    #initial stiffness of the p-y curve
    if phi < 30.0 /180*np.pi:
        kpy = 6.8  * 1.0e6/1.0e9 #in [N/mm3]
    elif phi < 36.0 /180*np.pi:
        kpy = 24.4 * 1.0e6/1.0e9
    else:
        kpy = 61.0 * 1.0e6/1.0e9    
    #equivalent depth
    Xreq2 = sigmav / p['gamma']
    #hyperbolic tangent curve parameters - chosen in such a way that the initial stiffness of the curve matched the suggested initial stiffness by Reeese&VanImple for medium dense dry sand
    #if min(kpy * Xreq2) < 0.001:
    #    print(min(kpy * Xreq2), kpy, min(sigmav), min(X))
    bl = pu * dsegmsingle['d'] / (kpy * Xreq2)  # bl    = (3.0*dsegmsingle['d']/80.0) / np.arctanh(0.95)   
    #axial mobilisation
    ba = 0.5 #90% of friction mobilised after 0.74mm displacement
    
    #mobiliation distances
    umob, wmob = func_mobilisationdistance(uext, X, Y, u, w, t, dsegmsingle, p)
    #mobilisation parameters
    zetaa, zetal = func_mobilisationparameters(umob, wmob, ba, bl)

    #AXIAL RESISTANCE
    #average normal stress
    sigman = sigmav*p['K'] + np.abs(zetal)*(0.5*sigmav*Ka - 0.5*sigmav*K0 + 0.25*np.abs(pu))
    #friction
    tau_u = sigman * np.tan(p['delta']*np.pi/180.0)

    #mobilsed lateral resistance    
    ql = dsegmsingle['d'] * pu * zetal
    #mobilsed axial resistance    
    qa = np.pi * dsegmsingle['d'] * tau_u * zetaa

    #return
    return(qa, ql)
    

#Function to calculate internal forces
def func_internalforcesnonincremental(dt, ddt, eps, dsegmsingle):
    """
    FUNCTION to calculate internal forces based on dimensional root state 
             forces are defined as parallel (N), perpendicular (V) to the DEFORMED root axis
    INPUT
    - <dt>:     dtheta/ds (dimensional, i.e. non-normalised)
    - <ddt>:    d^2theta/ds^2 (dimensional)
    - <eps>:    axial strain in s-direction (dimensional)
    - <dsegmsingle>: dictionary with segment data for one single segment
    OUTPUT
    - <out>:    dictionary with fields for axial force <N>, lateral force <N> and moment force <M>
    """
    #calculate strains
    eps_ax = eps
    eps_be = 0.5 * dt * dsegmsingle['d']
    #calculate secant stiffness
    Et = func_secantstiffness(eps_ax, dsegmsingle['Et1'], dsegmsingle['Et2'], dsegmsingle['Et3'])
    Eb = func_secantstiffness(eps_be, dsegmsingle['Eb1'], dsegmsingle['Eb2'], dsegmsingle['Eb3'])
    #Internal forces (axial N, shear V and moment M)
    N =  Et * dsegmsingle['A'] * eps
    V = -Eb * dsegmsingle['I'] * ddt
    M = -Eb * dsegmsingle['I'] * dt
    #return
    return(dict(zip(['N','V','M'], [N,V,M])))


#Function to calculate internal stresses
def func_internalstressnonincremental(dt, ddt, eps, dsegmsingle):
    """
    FUNCTION to calculate internal stresses and strains based on non-incremental, dimensional root state
    INPUT
    INPUT
    - <dt>:     dtheta/ds (dimensional, i.e. non-normalised)
    - <ddt>:    d^2theta/ds^2 (dimensional)
    - <eps>:    axial strain in s-direction (dimensional)
    - <dsegmsingle>: dictionary with segment data for one single segment
    OUTPUT
    - <out>:    dictionary with fields for axial force <N>, lateral force <N> and moment force <M>
    """
    #calculate strains
    eps_ax = eps
    eps_be = 0.5 * dt * dsegmsingle['d']
    #calculate secant stiffness
    Et = func_secantstiffness(eps_ax, dsegmsingle['Et1'], dsegmsingle['Et2'], dsegmsingle['Et3'])
    Eb = func_secantstiffness(eps_be, dsegmsingle['Eb1'], dsegmsingle['Eb2'], dsegmsingle['Eb3'])
    #calcaulte stresses
    sig_ax =  Et * eps
    sig_be = -Eb * (0.5*dsegmsingle['d'] * dt)
    sig_sh =  4.0/3.0 * (-Eb*dsegmsingle['I']/dsegmsingle['A']*ddt)
    return(dict(zip(['sig_ax', 'sig_be', 'sig_sh'], [sig_ax, sig_be, sig_sh])))


#Function to calculate internal stresses
def func_internalstrainnonincremental(dt, ddt, eps, dsegmsingle):
    """
    FUNCTION to calculate internal stresses and strains based on non-incremental, dimensional root state
    INPUT
    - <dt>:     dtheta/ds (dimensional, i.e. non-normalised)
    - <ddt>:    d^2theta/ds^2 (dimensional)
    - <eps>:    axial strain in s-direction (dimensional)
    - <dsegmsingle>: dictionary with segment data for one single segment
    OUTPUT
    - <out>:    dictionary with fields for axial force <N>, lateral force <N> and moment force <M>
    """
    eps_ax =  eps
    eps_be =  dt * 0.5*dsegmsingle['d']
    return(dict(zip(['eps_ax', 'eps_be'], [eps_ax, eps_be])))
        

#function to calculate global displacements
def func_globaldisplacements(u, w, t, Theta):
    """
    FUNCTION to calculate global displacements from local displacements
    INPUT
    - <u>:     local axial deformation
    - <w>:     local lateral deformation
    - <t>:     local rotation 
    - <Theta>: global orientaiton of non-displaced root segment
    OUTPUT
    - <out>:   dictionary with fields for global deformations <U>, <W> and <Theta>
    """
    U, W = mGeom.func_localdeformation2globaldeformations(u, w, Theta)
    return(dict(zip(['U', 'W', 'Theta'], [U, W, t])))


#function to calculate global displacements
def func_globalforces(N, V, M, Theta):
    """
    FUNCTION to calculate global internal forces from local internal forces
    INPUT
    - <N>:     axial force
    - <V>:     shear force
    - <M>:     moment
    - <Theta>: global orientaiton of non-displaced root segment
    OUTPUT
    - <out>:    dictionary with fields for global deformations <U>, <W> and <Theta>
    """
    FX, FY = mGeom.func_localdeformation2globaldeformations(N, V, Theta)
    return(dict(zip(['FX', 'FY', 'FTheta'], [FX, FY, M])))


def func_conditionatsegmentendsnonincremental(t, dt, ddt, eps, u, w, dsegmsingle, ends=True):
    """
    FUNCTION to calculate GLBOBAL, dimensional displacements and internal forces at both sides of a single root segment
    INPUT
    - <t>:         theta          (local, dimensional)
    - <dt>:        dtheta/ds      (local, dimensional)
    - <ddt>:       d^2theta/ds^2  (local, dimensional)
    - <eps>:       axial strain   (local, dimensional)
    - <u>:         u              (local, dimensional)    
    - <w>:         w              (local, dimensional)    
    - <dsegmsingle>: Dictionary wtih data for a single root segment
    OPTIONAL INPUT
    - <ends>         True  -> get forces acting on both ends of a segment (make forces on negative interface (left side) negative
                     False -> don't change orientation of forces
    OUTPUT
    - <out>:    Dictionary with output for U, W, Theta, FX, FY and FTheta at left and right side
                All are defined in the GLOBAL coorinate system, in dimenional terms
    """
    
    #output dictionary
    out = dict()
    
    #Displacements in GLOBAL direction (dimensional)
    out = func_globaldisplacements(u, w, t, dsegmsingle['Theta'])
    
    #Internal forces (axial N, shear V and moment M) - in direction of DEFORMED ROOT AXIS
    NVM = func_internalforcesnonincremental(dt, ddt, eps, dsegmsingle)
    #rotate to get forces in direction of non-deformed root
    Fx, Fy = mGeom.func_localdeformation2globaldeformations(NVM['N'], NVM['V'], t)
    Ftheta = NVM['M']
    #at first node, make internal forces negative since they work on the plane on the negative side of the segment
    if ends is True:
        Fx[0]     = -Fx[0]
        Fy[0]     = -Fy[0]
        Ftheta[0] = -Ftheta[0]

    #incremental Forces in GLOBAL direction
    out.update(func_globalforces(Fx, Fy, Ftheta, dsegmsingle['Theta']))
    #return
    return(out)    
    

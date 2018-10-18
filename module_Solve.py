# Solve - version 01/09/2017

#IMPORT PACKAGES
import numpy as np
import module_SoilRoot as mSoil
import module_Geometry as mGeom


# CREATE INITIAL STATE
def func_initialstate(dsegm, p, fields=['u','du','ddu','w','dw','ddw','dddw', 'us','ws']):
    """
    FUNCTION to create initial state of the root system (in local coordinate system, dimensionless parameters)
    INPUT
    - <dsegm>:     Root segment data
    - <p>:         Model input parameters
    OPTIONAL INPUT
    - <var_names>: List with variable names to be initialised, apart from <x>. All values are assumed 0.0
    OUTPUT
    - <v0>:        Dictionary with initial state of root
    """
    #initialise dictionary
    v0 = dict()
    #normalised results
    for i in dsegm.keys():
        v0[i] = dict(zip(fields, [np.zeros((1,p['n_node']))[0] for j in range(len(fields))]))
        v0[i]['x'] = np.linspace(0.0, 1.0, p['n_node'])
    #return
    return(v0)


# STEP SIZE
def func_stepsizeincrement(p):
    """
    FUNCTION to calculate soil increment size
    INPUT
    - <p>:         Model input parameters
    OUTPUT
    - <uext>:      Soil displacement increment
    """
    #Function to determine step size
    return(p['uext_inc0'])
    

#INITIAL GUESS FOR SOLVE_BVP SOLVER
def func_initialguess_initial(p):
    """
    FUNCTION to determine the initial guess for BVP_SOLVE (estimated in dimensionless coordinates).
             This function returns zero-vectors, so only usuable for the very first step
    INPUT
    - <p>:  Model input parameters
    OUTPUT
    - <init_guess_x, init_guess_y>: initial guess for x and y
    """
    #make guess
    init_guess_sh = np.linspace(0.0, 1.0, p['n_node'])
    init_guess_zh = np.zeros((p['n_var']*p['n_segm'], p['n_node']))
    #return guess
    return(init_guess_sh, init_guess_zh)


#INITIAL GUESS FOR SOLVE_BVP SOLVER
def func_initialguess_subsequent(sol_interp,  uext, uext_new):
    """
    FUNCTION to determine the initial guess for BVP_SOLVE (estimated in dimensionless coordinates).
             This function returns zero-vectors, so only usuable for the very first step
    INPUT
    - <p>:  Model input parameters
    OUTPUT
    - <init_guess_x, init_guess_y>: initial guess for x and y
    """
    zh_interp = np.zeros(np.shape(sol_interp[-1]))
    for i in range(zh_interp.shape[0]):
        for j in range(zh_interp.shape[1]):
            #zh_interp[i,j] = interpolate.CubicSpline(uext, [k[i,j] for k in sol_interp], bc_type="natural", extrapolate=True)(uext_new)
            ft = np.polyfit(uext, [k[i,j] for k in sol_interp], deg=1)
            zh_interp[i,j] = np.polyval(ft, uext_new)
    #np.polyfit(x, y, 3)
    return(zh_interp)
    

#SOLVE FUNCTION
def func_solve_wrapper(uext, dsegm, dnode, p):
    """
    FUNCTION to create function <func_solve> that calculates set of 1st order differentiation requried for bvp_solve
    INPUT
    - <uext>:   External displacement (dimensional coordinates)
    - <dsegm>:  Segment data
    - <dnode>:  Node data
    - <p>:      Model input parameters
    OUTPUT
    - <func_solve>: function that calculates first order derivates at location <s> (column position) for every variable (row) 
    """

    #solve function
    def func_solve(sh, zh):
        """
        FUNCTION that calculates the derivative of input vector <zh> at locations <sh>
		         all values are normalised by the segment length ('h' after parameter name indicates normalised parameter)
				 Thus all segments have the same normalised length, and the differential equaions for all segments can be solved simultaneously
        INPUT
        - <sh>: 1D array with non-dimensional positions along the displaced root axis (0 >= sh >= 1)
        - <zh>: 2D array with non-dimensional <z> at every location <x>
                zh = [th, dth/dsh, d^2th/dsh^2, epsh, uh, wh]  (all normalised/non-dimensional)
                when there are multiple segments, y is grouped as follows:
                zh = [th1, dth1/dsh, ..,  wh1,   th2, dth2/dsh ... , etc] etc
        - ...function uses many input parameters from functions with higher hierarchy
        OUTPUT
        - <dzh/dsh>: 1st order differentiation of <zh> with respect to <sh>
                     Size of the output array is similar to <zh.shape>
        """

        #function for single segment
        def func_solve_single(SegmentID):
            """
            FUNCTION to calculate first order differentiation for a single segment
            INPUT
            - <SegmentID>:  ID of segment for which calculation is to be performed
            - ...function uses many input parameters from functions with higher hierarchy
            OUTPUT
            - array with size <n_variables, n_node>, containing derivatives for all variables
            """            
            #Segment ID order number
            ind = list(dsegm.keys()).index(SegmentID)
            
            #solution values (input), make input more legible
            th    = zh[ind * p['n_var'] + 0, :] # theta(s)
            dth   = zh[ind * p['n_var'] + 1, :] # dtheta/ds 
            ddth  = zh[ind * p['n_var'] + 2, :] # d^2theta/ds^2
            epsh  = zh[ind * p['n_var'] + 3, :] # eps(s)
            uh    = zh[ind * p['n_var'] + 4, :] # u(s), normalised over length L, 
            wh    = zh[ind * p['n_var'] + 5, :] # w(s), normalised over length L, 
            
            #get LOCAL current root position - in non-normalised parameters
            xr = (sh + uh) * dsegm[SegmentID]['L']
            yr = (wh) * dsegm[SegmentID]['L']
            #get current GLOBAL current root position - in non-normalised parameters
            Xr, Yr = mGeom.func_localtoglobalsegmentpositions(xr, yr, dsegm[SegmentID]['X1'], dsegm[SegmentID]['Y1'], dsegm[SegmentID]['Theta'])

            #make parameters dimensional
            s = sh * dsegm[SegmentID]['L']
            u = uh * dsegm[SegmentID]['L']
            w = wh * dsegm[SegmentID]['L']
            t = th
    
            #external loading
            qa, ql = mSoil.func_soilresistance(uext, s, u, w, t, dsegm[SegmentID], p)
            
            #calculate strains - dimensional
            eps = epsh
            dt  = dth * dsegm[SegmentID]['L']**-1
            eps_ax = eps
            eps_be = 0.5 * dt * dsegm[SegmentID]['d']
            #calculate secant stiffness
            Et = mSoil.func_secantstiffness(eps_ax, dsegm[SegmentID]['Et1'], dsegm[SegmentID]['Et2'], dsegm[SegmentID]['Et3'])
            Eb = mSoil.func_secantstiffness(eps_be, dsegm[SegmentID]['Eb1'], dsegm[SegmentID]['Eb2'], dsegm[SegmentID]['Eb3'])
            #Segment stiffness parameters (normalised)
            EAL2 = Et * dsegm[SegmentID]['A'] / dsegm[SegmentID]['L']**2
            EIL4 = Eb * dsegm[SegmentID]['I'] / dsegm[SegmentID]['L']**4

            #divide soil resistances over segment length
            qlL = ql / dsegm[SegmentID]['L']
            qaL = qa / dsegm[SegmentID]['L']
           
            #calculate derivatives
            out_dth    = dth
            out_ddth   = ddth
            out_dddth  = 1.0/EIL4 * (qlL + EAL2*epsh*dth)
            out_depsh  = 1.0/EAL2 * (-EIL4*dth*ddth - qaL)
            out_duh    = (1.0 + epsh) * np.cos(th) - 1.0
            out_dwh    = (1.0 + epsh) * np.sin(th)
            
            #return array (rows=variables derivative, columns=values at x)
            return(np.vstack((out_dth, out_ddth, out_dddth, out_depsh, out_duh, out_dwh)))
        
        #stitch solutions for single segment together
        return(np.vstack(tuple([func_solve_single(ikey) for ikey in dsegm.keys()])))

    #return solver function    
    return(func_solve)


#FUNCTION TO CALCULATE ROOT STATE AT BOTH ENDS OF SEGMENT (incremental, displacements and internal forces in global dimenional coordinate system
def func_bvpconditionsatrootends(zha, zhb, dsegm, p):
    """
    FUNCTION to calculate GLOBAL, dimensional displacements and internal forces at both sides of all root segments, based on <solve_bvp> boundary format
    INPUT
    - <zha>:    solve_bvp solution at left side of segment  (NodeID1) (dimensionless)
    - <zhb>:    solve_bvp solution at right side of segment (NodeID2) (dimensionless)
    - <dsegm>:  Segment data
    - <p>:      Model input parameters
    OUTPUT
    - <out>:    Dictionary with, per segment, output for U, W, Theta, FX, FY and FTheta at left and right side
                All are defined in the GLOBAL coorinate system, in dimenional terms
    """
       
    #combine ya and yb into single array
    zh = np.stack((zha, zhb), axis=1)
    
    #get values for every parameter from output
    th    = zh[np.arange(0, p['n_segm'])*p['n_var']+0, :] #zh[[i*p['n_var']+0 for i in range(p['n_segm'])], :]
    dth   = zh[np.arange(0, p['n_segm'])*p['n_var']+1, :] #zh[[i*p['n_var']+1 for i in range(p['n_segm'])], :]
    ddth  = zh[np.arange(0, p['n_segm'])*p['n_var']+2, :] #zh[[i*p['n_var']+2 for i in range(p['n_segm'])], :]
    epsh  = zh[np.arange(0, p['n_segm'])*p['n_var']+3, :] #zh[[i*p['n_var']+3 for i in range(p['n_segm'])], :]
    uh    = zh[np.arange(0, p['n_segm'])*p['n_var']+4, :] #zh[[i*p['n_var']+5 for i in range(p['n_segm'])], :]
    wh    = zh[np.arange(0, p['n_segm'])*p['n_var']+5, :] #zh[[i*p['n_var']+6 for i in range(p['n_segm'])], :]

    #segment lengths
    L   = np.array([i['L'] for sID,i in dsegm.items()])
    #get normalised values
    t   = th   * L[:, np.newaxis]**0
    dt  = dth  * L[:, np.newaxis]**-1
    ddt = ddth * L[:, np.newaxis]**-2
    eps = epsh * L[:, np.newaxis]**0
    u   = uh   * L[:, np.newaxis]**1
    w   = wh   * L[:, np.newaxis]**1
        
    #output dictionary
    out = dict()
    
    #loop through segments
    for sID,j in zip(dsegm.keys(), range(p['n_segm'])):
        
        #create empty dictionary for each segment in output dictionary
        out[sID] = dict()
        
        #call function
        out[sID] = mSoil.func_conditionatsegmentendsnonincremental(t[j,:], dt[j,:], ddt[j,:], eps[j,:], u[j,:], w[j,:], dsegm[sID], ends=True)
        
    #return
    return(out)


#BOUNDARY CONDITION FUNCTION
def func_boundary_wrapper(uext, dsegm, dnode, p):
    """
    FUNCTION to create function <func_solve> that calculates set of 1st order differentiation requried for bvp_solve
    INPUT
    - <uext>:      External displacement (dimensional coordinates)
    - <dsegm>:     Segment data
    - <dnode>:     Node data
    - <p>:         Model input parameters
    OUTPUT
    - <func_boundary>: function that calculates boundary conditions based on bvp_solution at left and right side. 
                       target value of all bc's is 0.0
    """
    
    def func_boundary(zha, zhb):
        """
        FUNCTION that defines boundary conditions
        INPUT
        - <zha>: solve_bvp solution at left boundary,  normalised parameters
        - <zhb>: solve_bvp solution at right boundary, normalised parameters
        OUTPUT
        - <bc>: list with values that should be all equal to 0.0
                Length of the 1D output array is similar to <zha.shape>
        """
              
        #initiate bc output list
        bc = list()  
        
        #conditions at ends of root segments
        cend = func_bvpconditionsatrootends(zha, zhb, dsegm, p)
        #soil deformations at original positions of node - assumes nodes are fixed in X and Y directions
        NodePos = mGeom.func_nodaldeformations(uext, dnode, p)
        
        #loop through all nodes        
        for nID,i in dnode.items():
            
            #X-direction
            if i['bound_X'] == True:
                #clamped end, no displacement
                for sID,sSide in zip(i['SegmentID'], i['SegmentSide']):
                    bc.append(cend[sID]['U'][sSide] - NodePos[nID]['Us'])
            else:
                #free end, residual force is zero
                bc.append(sum([cend[sID]['FX'][sSide] for sID,sSide in zip(i['SegmentID'], i['SegmentSide'])]))
                #free end, displacements of connecting nodes should all be equal
                if i['NodeType'] == 'middle':
                    #middle node, u1=u2, u2=u3 etc
                    for sID1,sSide1,sID2,sSide2 in zip(i['SegmentID'][0:-1], i['SegmentSide'][0:-1], i['SegmentID'][1:], i['SegmentSide'][1:]):
                        bc.append(cend[sID1]['U'][sSide1] - cend[sID2]['U'][sSide2])

            #Y-direction
            if i['bound_Y'] == True:
                #clamped end, no displacement
                for sID,sSide in zip(i['SegmentID'], i['SegmentSide']):
                    bc.append(cend[sID]['W'][sSide] - NodePos[nID]['Ws'])
            else:
                #free end, residual force is zero
                bc.append(sum([cend[sID]['FY'][sSide] for sID,sSide in zip(i['SegmentID'], i['SegmentSide'])]))
                #free end, displacements of connecting nodes should all be equal
                if i['NodeType'] == 'middle':
                    #middle node, u1=u2, u2=u3 etc
                    for sID1,sSide1,sID2,sSide2 in zip(i['SegmentID'][0:-1], i['SegmentSide'][0:-1], i['SegmentID'][1:], i['SegmentSide'][1:]):
                        bc.append(cend[sID1]['W'][sSide1] - cend[sID2]['W'][sSide2])
            
            #Theta-direction
            if i['bound_Theta'] == True:
                #clamped end, no displacement
                for sID,sSide in zip(i['SegmentID'], i['SegmentSide']):
                    bc.append(cend[sID]['Theta'][sSide])
            else:
                #free end, residual force is zero
                bc.append(sum([cend[sID]['FTheta'][sSide] for sID,sSide in zip(i['SegmentID'], i['SegmentSide'])]))
                #free end, displacements of connecting nodes should all be equal
                if i['NodeType'] == 'middle':
                    #middle node, u1=u2, u2=u3 etc
                    for sID1,sSide1,sID2,sSide2 in zip(i['SegmentID'][0:-1], i['SegmentSide'][0:-1], i['SegmentID'][1:], i['SegmentSide'][1:]):
                        bc.append(cend[sID1]['Theta'][sSide1] - cend[sID2]['Theta'][sSide2])
        
        #return array with boundary conditions - no deform left - free end right
        return(np.array(bc)) 
        
    #return wrapper output - function        
    return(func_boundary)
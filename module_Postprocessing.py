# Postprocessing functions - version 01/09/2017

#IMPORT PACKAGES
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import cm
import module_Geometry as mGeom
import module_SoilRoot as mSoil


#CONVERT BVP_SOLVE RESULTS TO DICTIONARY FOR INCREASED READABILITY
def func_sol2dict(sol, dsegm, p):
    """
    FUNCTION to convert bvp_solve function results to dictionary - dimensionless results
    INPUT
    - <soil>:  Solution array, size = (n_segm*n_var) x (n_node)
    - <dsegm>: Segment data
    - <p>:     Model input parameters
    OUTPUT
    - <out>:   dictionary with solve_bvp results per segment per variable
    """
    #create output dictionary
    out = dict()
    #loop through segments
    for sID,i in zip(dsegm.keys(), range(p['n_segm'])):
        out[sID] = dict()
        out[sID]['sh']   = np.linspace(0.0, 1.0, p['n_node']) #Normalised coordinate along the root axis
        out[sID]['th']   = sol[i*p['n_var']+0, :]
        out[sID]['dth']  = sol[i*p['n_var']+1, :]
        out[sID]['ddth'] = sol[i*p['n_var']+2, :]
        out[sID]['epsh'] = sol[i*p['n_var']+3, :]
        out[sID]['uh']   = sol[i*p['n_var']+4, :]
        out[sID]['wh']   = sol[i*p['n_var']+5, :]
        out[sID]['xh']   = out[sID]['sh'] + out[sID]['uh']
        out[sID]['yh']   = out[sID]['wh']
    #return
    return(out)
 
    
#CONVERT SOLVE RESULTS FROM DIMENSIONLESS INTO DIMENSIONAL RESULTS - USE OUTPUT FROM <func_sol2dict>
def func_sol2soldimensional(sol, dsegm, p):
    """
    FUNCTION to convert bvp_solve function results to dictionary - dimensional results
    INPUT
    - <soil>:  Solution array, size = (n_segm*n_var) x (n_node)
    - <dsegm>: Segment data
    - <p>:     Model input parameters
    OUTPUT
    - <out>:   dictionary with solve_bvp results per segment per variable
    """
    #create output dictionary
    out = dict()
    #loop through segments
    for sID,i in zip(dsegm.keys(), range(p['n_segm'])):
        out[sID] = dict()
        out[sID]['s']   = np.linspace(0.0, dsegm[sID]['L'], p['n_node']) #Normalised coordinate along the root axis
        out[sID]['t']   = sol[i*p['n_var']+0, :] * dsegm[sID]['L']**0
        out[sID]['dt']  = sol[i*p['n_var']+1, :] * dsegm[sID]['L']**-1
        out[sID]['ddt'] = sol[i*p['n_var']+2, :] * dsegm[sID]['L']**-2
        out[sID]['eps'] = sol[i*p['n_var']+3, :] * dsegm[sID]['L']**0
        out[sID]['u']   = sol[i*p['n_var']+4, :] * dsegm[sID]['L']**1
        out[sID]['w']   = sol[i*p['n_var']+5, :] * dsegm[sID]['L']**1
        out[sID]['x']   = out[sID]['s'] + out[sID]['u']
        out[sID]['y']   = out[sID]['w']
        out[sID]['Xr'], out[sID]['Yr'] = mGeom.func_localtoglobalsegmentpositions(out[sID]['x'], out[sID]['y'], dsegm[sID]['X1'], dsegm[sID]['Y1'], dsegm[sID]['Theta'])
        
    #return
    return(out)    


#COMBINE DICTIONARIES CONTAINING SEGMENT RESULTS
def func_combineresults(d1, d2):
    """
    FUNCTION to combine dictionaries containing segmental data/results. 
    INPUT
    - <d1,d2>:  Dictionary with segmentID's as keys. For every segment, it contains another dictionary with data, using variable names as keys
    OUTPUT
    - <out>:    Combined dictionary
    """
    #output dictionary
    out = d1.copy()
    #loop through segments
    for sID in d2.keys():
        out[sID].update(d2[sID])
    #return
    return(out)
    

#CALCULATE INTERNAL FORCES, STRAINS AND STRESSES (NON_INCREMENTAL)
def func_internalforces(v, dsegm):
    """
    FUNCTION to calculate internal root forces, strains and stresses (dimensional, non-incremental)
    INPUT
    - <v>:      Current root state (dimensional, local coordinate system)
    - <dsegm>:  Segment data
    OUTPUT
    - <out>:    Combined dictionary
    """
    #output dictionary
    out = dict()
    #loop through segments
    for sID,i in v.items():
        #emtpy dictionary
        out[sID] = dict()
        #Internal root strains
        out[sID] =      mSoil.func_internalstrainnonincremental(i['dt'], i['ddt'], i['eps'], dsegm[sID])
        #Internal root forces
        out[sID].update(mSoil.func_internalforcesnonincremental(i['dt'], i['ddt'], i['eps'], dsegm[sID]))
        #Internal root stresses
        out[sID].update(mSoil.func_internalstressnonincremental(i['dt'], i['ddt'], i['eps'], dsegm[sID]))
    #return
    return(out)


#CALCULATE REINFORCEMENT PER SEGMENT
def func_reinforcementspersegment(F, v, dsegm, p):
    """
    FUNCTION to calcualte reinforcing forces in global X and Y direction, for every point in every segment
    INPUT
    - <F>:      Current internal forces and stresses (dimensional, aligned with displaced root axis)
    - <v>:      Current root state (dimensional, local)
    - <dsegm>:  Segment data
    - <p>:      Model input parametrs
    OUTPUT
    - <out>:    Dictionary with reinforcements per segment 
                FX = Force in X-direction
                FY = Force in Y-direction
                F  = Total force, F=Fy+Fx*tan(critical state friction angle)
    """
    #output dictionary
    out = dict()
    #loop through segments
    for sID,i in F.items():
        #emtpy dictionary
        out[sID] = dict()
        #Get internal forces in global coordinate system
        Fx, Fy = mGeom.func_localdeformation2globaldeformations(i['N'], i['V'], dsegm[sID]['Theta']+v[sID]['t'])
        #Reinforcing forces in global X and Y directions
        out[sID]['FX'] = Fx
        out[sID]['FY'] = Fy
        #total reinforcement
        out[sID]['F']  = out[sID]['FY'] + out[sID]['FX'] * np.tan(p['phi_cs']/180.0*np.pi)
    return(out)


#FUNCTION TO CALCULATE TOTAL REINFORCEMENT AT EACH DEPTH LEVEL
def func_reinforcementsperdepth(V, F, n=100):
    """
    FUNCTION to calcualte reinforcing forces in global X and Y direction, for every point in every segment
    INPUT
    - <V>:      Current root state (dimensional, global)
    - <F>:      Current reinforcements per segment
    OPTIONAL INPUT
    - <n>:      Number of points, linearly distributed over range of X
    OUTPUT
    - <out>:    Dictionary with reinforcements per segment 
                FX = Force in X-direction
                FX = Force in Y-direction
                F  = Total force, F=FY+FX*tan(critical state friction angle)
    """
    #output dictionary
    out = dict()
    #Get X-corodinates at which to interpolate results
    out['X'] = np.linspace(np.min([i['Xr'] for sID,i in V.items()]), np.max([i['Xr'] for sID,i in V.items()]), n)
    #loop through segments
    for sID,i in V.items():
        if sID == list(V.keys())[0]:
            out['FX']  = interpolate.interp1d(i['Xr'], F[sID]['FX'], kind='linear', bounds_error=False, fill_value=(0,0))(out['X'])
            out['FY']  = interpolate.interp1d(i['Xr'], F[sID]['FY'], kind='linear', bounds_error=False, fill_value=(0,0))(out['X'])
            out['F']   = interpolate.interp1d(i['Xr'], F[sID]['F'] , kind='linear', bounds_error=False, fill_value=(0,0))(out['X'])
        else:
            out['FX'] += interpolate.interp1d(i['Xr'], F[sID]['FX'], kind='linear', bounds_error=False, fill_value=(0,0))(out['X'])
            out['FY'] += interpolate.interp1d(i['Xr'], F[sID]['FY'], kind='linear', bounds_error=False, fill_value=(0,0))(out['X'])
            out['F']  += interpolate.interp1d(i['Xr'], F[sID]['F'] , kind='linear', bounds_error=False, fill_value=(0,0))(out['X'])  
    return(out) 


def func_reinforcementatshearplane(Fr, v, dsegm, p):
    """
    FUNCTION to calcaulte reinforcements at shear plane, for all segments individually
    INPUT
    - <Fr>:     dictionary with reinforcements per segment (dimensional, in GLOBAL coordinates)
    - <v>:      dictionary containing root postions (global, dimensional), in fields <Xr> and <Yr>
    - <dsegm>:  segment data
    - <p>:      model input parameters
    OUTPUT
    - <out>:    Dictionary with reinforcements per segment 
                Fdirect = Force in X-direction
                Fconfine = Force in Y-direction
                F  = Total force, F=FY+FX*tan(critical state friction angle)
    """
    #output dictionary
    out = dict()
    #loop through segments
    for sID,i in Fr.items():
        #root positions in shear plane coordinate system (x=aligned with shear plane orientation, y=perpendicular)
        [xsh, ysh] = mGeom.func_globaldeformation2localdeformations(v[sID]['Xr']-p['shearplane_depth'], v[sID]['Yr'], p['shearplane_angle'])
        #linear interpolation of global forces at intersection
        FXsh = interpolate.interp1d(ysh, i['FX'], kind='linear', bounds_error=False, fill_value=0.0)(0.0)
        FYsh = interpolate.interp1d(ysh, i['FY'], kind='linear', bounds_error=False, fill_value=0.0)(0.0)
        #express in orientation of shear plane
        [F_confine, F_direct] = mGeom.func_globaldeformation2localdeformations(FXsh, FYsh, np.pi/2.0-p['shearplane_angle'])
        out[sID] = dict(zip(['F_direct', 'F_confine', 'F_total'], [F_direct, F_confine, F_direct+F_confine*np.tan(p['phi_cs']/180.0*np.pi)]))
    return(out)
        

def func_maxstresses(v, F, dsegm, p):
    """
    FUNCTION to calcaulte maximum stresses in segments 
    INPUT
    - <v>:      displaced root state
    - <F>:      dictionary with internal forces and stresses per segment (dimensional, in coordinates local to deformed root axis)
    - <dsegm>:  segment data
    - <p>:      model input parameters
    OUTPUT
    - <out>:    Dictionary with max stresses
                <sig_be_shearplane> = maximum bending stress at shear plane
                <sig_ax_shearplane> = maximum tensile stress at shear plane
                <sig_sh_shearplane> = maximum shear   stress at shear plane
    """
    #output dictionary
    out = dict()
    #loop through segments
    for sID,i in F.items():
        #create dictionary for each segment
        out[sID] = dict()
        #root positions in shear plane coordinate system (x=aligned with shear plane orientation, y=perpendicular)
        [xsh, ysh] = mGeom.func_globaldeformation2localdeformations(v[sID]['Xr']-p['shearplane_depth'], v[sID]['Yr'], p['shearplane_angle'])
        #stresses at shear plane
        out[sID]['sig_ax_shearplane'] = np.asscalar(interpolate.interp1d(ysh, i['sig_ax'], kind='linear', bounds_error=False, fill_value=0.0)(0.0))
        out[sID]['sig_be_shearplane'] = np.asscalar(interpolate.interp1d(ysh, i['sig_be'], kind='linear', bounds_error=False, fill_value=0.0)(0.0))
        out[sID]['sig_sh_shearplane'] = np.asscalar(interpolate.interp1d(ysh, i['sig_sh'], kind='linear', bounds_error=False, fill_value=0.0)(0.0))
        #maximum stresses
        out[sID]['sig_ax_max'] = np.max(i['sig_ax'])
        out[sID]['sig_be_max'] = np.max(i['sig_be'])
        out[sID]['sig_sh_max'] = np.max(i['sig_sh'])
        #location of maximum stresses
        out[sID]['X_sig_ax_max'] = v[sID]['Xr'][np.argmax(i['sig_ax'])]
        out[sID]['X_sig_be_max'] = v[sID]['Xr'][np.argmax(i['sig_be'])]
        out[sID]['X_sig_sh_max'] = v[sID]['Xr'][np.argmax(i['sig_sh'])]
        #maximum fibre stress
        sig_fi = i['sig_ax'] + np.abs(i['sig_be'])
        out[sID]['sig_fi_max']    = np.max(sig_fi)
        out[sID]['sig_fi_max_ax'] =     F[sID]['sig_ax'][np.argmax(sig_fi)]
        out[sID]['sig_fi_max_be'] = abs(F[sID]['sig_be'][np.argmax(sig_fi)])
        out[sID]['X_sig_fi_max']  = v[sID]['Xr'][np.argmax(sig_fi)]
    return(out)


def mobilisationdistances(v, uext, dsegm, p):
    for sID,i in v.items():
        i['umob'], i['wmob'] = mSoil.func_mobilisationdistance(uext, i['Xr'], i['Yr'], i['u'], i['w'], i['t'], dsegm[sID], p)


def func_plotsingle(v, Vr, Vs, F, dsegm, p, size=(15.0, 9.0)):
    """ 
    FUNCTION that plots the results of a single calculation step
    INPUT
    - <v>:  dictionary with current root state (local coordinate system)
    - <Vr>: dictionary with current root position (global coordinate system)
    - <Vs>: dictionary with current soil position (global coordinate system)
    - <F>:  dictionary with current internal forces and stresses (coordinate system local to displaced root axis)
    - <p>:  dictionary with input parameters
    - <dsegm>: segment data
    OPTIONAL INPUT
    - <size> tuple with figure size (inches)
    OUTPUT
    - <f>: figure handle
    - <sp_array> figure subplot handles
    """

    #number of subplot rows and columns
    n_sp  = 9
    n_row = 3
    n_col = 3
    plt_count = -1
        
    #create subplots
    f, sp_array = plt.subplots(n_row, n_col, figsize=size)

    #subplot positions
    sp_pos = [(int(np.floor(i/n_col)), i%n_col) for i in range(n_sp)]

    #Segment colors
    colo_root = dict(zip(dsegm.keys(), [cm.jet(i) for i in np.linspace(0.0, 1.0, p['n_segm'])]))
    colo_root_orig = '0.5'

    #local, dimensional undisplaced root positions
    x = np.linspace(0, 1, p['n_node'])

    #PLOT - SHAPES
    plt_count = plt_count + 1
    for sID,i in Vr.items():
        #displaced root
        sp_array[sp_pos[plt_count]].plot(i['Xr'], i['Yr'],   
            color=colo_root[sID], 
            linestyle='-')
        #displaced soil
    for sID,i in Vs.items():    
        sp_array[sp_pos[plt_count]].plot(i['Xs'], i['Ys'],   
            color=colo_root[sID], 
            linestyle=':',
            linewidth=0.5)      
    #labels
    sp_array[sp_pos[plt_count]].set_xlabel('$X$ [mm]')
    sp_array[sp_pos[plt_count]].set_ylabel('$Y$ [mm]')

    #PLOT - AXIAL DISPLACEMENT
    plt_count = plt_count + 1
    for sID,i in v.items():
        sp_array[sp_pos[plt_count]].plot(x, i['u'],
            color=colo_root[sID], 
            linestyle='-')
    #labels
    sp_array[sp_pos[plt_count]].set_xlabel('x [-]')
    sp_array[sp_pos[plt_count]].set_ylabel('Axial displacement [mm]')

    #PLOT - LATERAL DISPLACEMENT
    plt_count = plt_count + 1
    for sID,i in v.items():
        sp_array[sp_pos[plt_count]].plot(x, i['w'],
            color=colo_root[sID], 
            linestyle='-')
    #labels
    sp_array[sp_pos[plt_count]].set_xlabel('x [-]')
    sp_array[sp_pos[plt_count]].set_ylabel('Lateral displacement [mm]')

    #PLOT - AXIAL STRESSES
    plt_count = plt_count + 1
    for sID,i in F.items():
        sp_array[sp_pos[plt_count]].plot(x, i['sig_ax'], 
            color=colo_root[sID], 
            linestyle='-')
    #labels
    sp_array[sp_pos[plt_count]].set_xlabel('x [-]')
    sp_array[sp_pos[plt_count]].set_ylabel('Axial stress [MPa]')

    #PLOT - BENDING STRESSES
    plt_count = plt_count + 1
    for sID,i in F.items():
        sp_array[sp_pos[plt_count]].plot(x, i['sig_be'],
            color=colo_root[sID], 
            linestyle='-')
    #labels
    sp_array[sp_pos[plt_count]].set_xlabel('x [-]')
    sp_array[sp_pos[plt_count]].set_ylabel('Bending stress [MPa]')
    
    #PLOT - SHEAR STRESSES
    plt_count = plt_count + 1
    for sID,i in F.items():
        sp_array[sp_pos[plt_count]].plot(x, i['sig_sh'], 
            color=colo_root[sID], 
            linestyle='-')
    #labels
    sp_array[sp_pos[plt_count]].set_xlabel('x [-]')
    sp_array[sp_pos[plt_count]].set_ylabel('Shear stress [MPa]')

    #PLOT - AXIAL MOBILISATION
    plt_count = plt_count + 1
    for sID,i in v.items():
        sp_array[sp_pos[plt_count]].plot(x, i['umob'], 
            color=colo_root[sID], 
            linestyle='-')
    #labels
    sp_array[sp_pos[plt_count]].set_xlabel('x [-]')
    sp_array[sp_pos[plt_count]].set_ylabel('Axial mobilisation [mm]')
    
    #PLOT - LATERAL MOBILISATION
    plt_count = plt_count + 1
    for sID,i in v.items():
        sp_array[sp_pos[plt_count]].plot(x, i['wmob'], 
            color=colo_root[sID], 
            linestyle='-')
    #labels
    sp_array[sp_pos[plt_count]].set_xlabel('x [-]')
    sp_array[sp_pos[plt_count]].set_ylabel('Lateral mobilisation [mm]')
    #right layout
    plt.tight_layout()

    #PLOT - ANGLES
    plt_count = plt_count + 1
    for sID,i in v.items():
        sp_array[sp_pos[plt_count]].plot(x, i['t'], 
            color=colo_root[sID], 
            linestyle='-')
    #labels
    sp_array[sp_pos[plt_count]].set_xlabel('x [-]')
    sp_array[sp_pos[plt_count]].set_ylabel('Theta [rad]')
    #right layout
    plt.tight_layout()

    #return
    return(f, sp_array)
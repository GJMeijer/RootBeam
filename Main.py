#version 20170831
#(c) G.J. Meijer, 2017

## USER INPUT SETTINGS

#file name with run parameters
file_name_parameters = 'Input_Parameters.csv'
#settings for plotting data, saving data etc
data_save = True
plot_show = False
plot_save = False


## MAIN CODE

#load packages
import os
from scipy.integrate import solve_bvp
import module_InputOutput as mIO
import module_SoilRoot as mSoil
import module_Solve as mSolve
import module_Postprocessing as mPost
import module_Geometry as mGeom

#change work directory to current directory
abspath = os.path.abspath(__file__)
directory = os.path.dirname(abspath)
os.chdir(directory)

#Load input parameters
print('LOAD INPUT PARAMETERS')
par = mIO.loadcsv(file_name_parameters, dict1=True, dict2=True)

#loop through input parameters - runs
for pkey,p in par.items():

    #print progress
    print('RUN ' + str(pkey) + ' out of ' + str(len(par)))
    
    #Load node and segment data    
    dnode, dsegm = mIO.loadnodesegment(p['NodeFile'], p['SegmentFile'])
    #add number/counts to p
    p['n_segm']     = len(dsegm)             #number of root segments
    p['n_node']     = int(p['n_node'])       #number of nodes per segment, convert to integer
    p['n_step_max'] = int(p['n_step_max'])   #max number of displaecment steps, convert to integer
    p['n_var']      = 6 					 #number of variables to solve per segment (th,dth,ddth,epsh,uh,wh)-->(dth,ddth,dddth,depsh,duh,dwh)

    #current displacement - initialise list
    uext = [0.0]
    #initialise shear reinforcement results dictionary
    Fsh = dict()
    #initialise stress results dictionary
    Fsig = dict()
    
    #initial step counter
    step_counter = 1
    #initial step size
    uext_increment = p['uext_inc0']
    #switch to indicate whether a next step should be calculated
    step_continue_switch = True
    #keep track of numbers of iterations
    n_iter_record = list()
    
    #loop through displacement steps
    #Continue looping when
    # - max step number not reached
    # - max shear deformation not reached
    # - boolean switch <step_continue_switch> is still true
    while step_counter <= p['n_step_max'] and uext[-1] <= p['uext_max'] and step_continue_switch == True: 

        #if first step - initialise current root state <v> and soil resistance state <q>
        if step_counter==1:
            #initial root position (global coordinate system, dimensional)
            V0 = mGeom.func_initialglobalrootposition(dsegm, p)
            #initial guess - normalised parameters (normalised by segment length)
            init_guess_sh, init_guess_zh = mSolve.func_initialguess_initial(p)
        
        #iteration parameters
        iteration_solved_switch = False
        iter_counter        = 1
        
        #loop until a suitable step size has been found
        while iteration_solved_switch is False and iter_counter <= p['n_iter_max']:
   
            #get step size
            if iter_counter == 1 :
                if step_counter >= 5 and set(n_iter_record[-5:])=={1}:
                    #increment can be increased - because the same for last 2 iteration and solution stable
                    uext_increment = p['uext_factor_increase'] * uext_increment                
            else:
				#displcement step too large - no solution found - reduce step size
                uext_increment = p['uext_factor_decrease'] * uext_increment
            #check if step size larger than minimum step size specified
            if uext_increment <= p['uext_incmin']:
				#Required step size is below user-defined minimum step size --> stop analysis
                iteration_solved_switch = True
                step_continue_switch = False
            else:
                #Step size acceptable, but cap with maximum user-defined step size
                uext_increment = min([uext_increment, p['uext_incmax']])
                #get current displacment <uext_temp>
                uext_temp = uext[-1] + uext_increment
                            
            #print progress
            print('- Step ' + str(step_counter) + '/' + str(p['n_step_max']) + ', uext ' + str(round(uext_temp,3)) + '/' + str(p['uext_max']) + ', iter ' + str(iter_counter) + '/' + str(int(p['n_iter_max'])))
                
            #Create functions for differential equation solver to solve, based on current parameters
    			#create solve function
            f_solve = mSolve.func_solve_wrapper   (uext_temp, dsegm, dnode, p)
            #create function to account for boundary conditions
            f_bound = mSolve.func_boundary_wrapper(uext_temp, dsegm, dnode, p)
    
            #solve set of differential equations
            sol_temp = solve_bvp(f_solve, f_bound, init_guess_sh, init_guess_zh, p=None, S=None, verbose=0, tol=p['solve_tolerance'])
           
            #decide what to do based on whether solver converges
            if sol_temp.success is True:
                #solution has been found - solver converges
                
                #set solved switch to true
                iteration_solved_switch = True
                
                #add current displacement to list of displacements
                uext.append(uext_temp)                
                
                #interpolate solution at certain points along the segment
                sol = sol_temp.sol(init_guess_sh)
                
                #transform results in more readable dictionary form - LOCAL coordinates, dimensional
                v = mPost.func_sol2soldimensional(sol, dsegm, p)
                
                #get internal forces
                F  = mPost.func_internalforces(v, dsegm)
                #get current root position in GLOBAL, dimensional coordinates
                Vr = mGeom.func_localdimensionless2globalpositions(v, dsegm)
                #get current soil position in GLOBAL, dimensional coordinates
                Vs = mSoil.func_currentglobalsoilposition(V0, uext_temp, p)
                
                #Get reinforcements per segment
                Fr = mPost.func_reinforcementspersegment(F, v, dsegm, p) 
                
                #add mobiilsation distances
                mPost.mobilisationdistances(v, uext_temp, dsegm, p)
                
                #Values on shear plane
                Fsh[step_counter] = mPost.func_reinforcementatshearplane(Fr, v, dsegm, p)   
                
                #Max stresses in segment
                Fsig[step_counter] = mPost.func_maxstresses(v, F, dsegm, p)
                
                #Create new initial guess based on last solution
                init_guess_zh = sol
                
                #Save step data
                if data_save is True and (step_counter-1)%p['SaveStep']==0:
                    #name of directory to save in
                    directory = 'Results/' + p['ModelName'] + '/Run' + str(pkey) + '/Root'
                    #save all step output data - root
                    file_name1 = p['ModelName'] + '_Results_RootStep' + str(step_counter) + '.csv'
                    mIO.dict2csv(file_name1, mIO.dictdict2dict(v, key1_header='SegmentID'), directory=directory)
        
                    #name of directory to save in
                    directory = 'Results/' + p['ModelName'] + '/Run' + str(pkey) + '/Soil'
                    #save all step output data - soil
                    file_name1 = p['ModelName'] + '_Results_SoilStep' + str(step_counter) + '.csv'
                    mIO.dict2csv(file_name1, mIO.dictdict2dict(Vs, key1_header='SegmentID'), directory=directory)
                
                #record number of iterations
                n_iter_record.append(iter_counter)
                
                #increase step counter
                step_counter += 1
                
            else:
                #solution has not been found
                
                #update iteration counter
                iter_counter += 1
                
                #kill analysis if last iteration - break out of step loop because no suitable solution found to continue with
                if iter_counter == p['n_iter_max']:
                    #break iteration loop
                    iteration_solved_switch = True
                    #break step loop
                    step_continue_switch = False
                
    #after last step
    if data_save is True:
        #save data
		
        #name of save directory
        directory = 'Results/' + p['ModelName'] + '/Run' + str(pkey) + '/Summary'
    
        #save all step output data - contributions to reinforcement
        file_name2 = p['ModelName'] + '_Results_ShearReinforcement.csv'
        mIO.dict2csv(file_name2, mIO.dictdict2dict({ikey:mIO.dictdict2dict(i, key1_header='SegmentID') for ikey,i in Fsh.items()}, key1_header='StepID'), directory=directory)
             
        #save force-displacement data
        file_name3 = p['ModelName'] + '_Results_ShearForceDisp.csv'
        Ftot = [0.0] + [sum([j['F_total'] for jkey,j in i.items()]) for ikey,i in Fsh.items()]
        isteps = list(range(0,step_counter))
        dic  = dict(zip(['StepID', 'uext', 'F_total'], [isteps, uext, Ftot]))
        mIO.dict2csv(file_name3, dic, directory=directory)
        
        #save data for stresses
        file_name4 = p['ModelName'] + '_Results_Stresses.csv'
        mIO.dict2csv(file_name4, mIO.dictdict2dict({ikey:mIO.dictdict2dict(i, key1_header='SegmentID') for ikey,i in Fsig.items()}, key1_header='StepID'), directory=directory)
        
    #plot final position
    if plot_show is True or plot_save is True:
		#create plot with final root position
        plot_handle, sp_array = mPost.func_plotsingle(v, Vr, Vs, F, dsegm, p, size=(9.0, 7.0))
        #save plot
        if plot_save == True:        
            directory = 'Results/' + p['ModelName'] + '/Run' + str(pkey) + '/Summary'
            plot_handle.savefig(directory + '/' + p['ModelName']+'_run'+str(pkey)+'_plotlaststep.pdf')
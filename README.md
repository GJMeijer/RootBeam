# RootBeam
Analysis of coupled axial and lateral deformation of roots in soil

Notes 27/09/2018, GJMeijer  
Python code to analyse the behaviour of roots loaded in direct shear by solving the (large deformation) bending equation, taking i to account both axial and bending stiffness.

## Publication details
This code belong to journal publication:  

> Analysis of coupled axial and lateral deformation of roots in soil  
> GJ Meijer, D Muir Wood, JA Knappett, AG Bengough, T Liang  
> International Journal for Analytical and Numerical Methods in Geomechanics  
> 2018, volume x, issue x, page xx-xx  
> DOI: xxx

This code can be used to reproduce some of the shear tests with ABS described in this paper (Liang et al.)


## Quickstart
Ensure that all python files (`.py`) and input files (`.csv`) are placed in the same directory.
Analysis can be run by executing the `Main.py` script.
In the first lines of `Main.py`, some user-defined settings are required
* The name of the input parameter file, including extension
* Some booleans to decide whether results should be plotted, saved etc.

## Python requirements
Code is written in Python 3
Executing required installed modules (apart from standard modules):
* `SciPy`
* `NumPy`
* `MatPlotLib`


# Notes on data format

## Units
Throughout the code, units used:
* Length in millimeters (mm)
* Force in newtons (N)
* Angles in radians (rad)
Note that therefore, in constrast to common geotech practice:
* Stesses and stiffnesses are in (MPa)
* Unit weights are in (N/mm3)

## Input parameter file - required parameters
Required columns:

| Parameter | Type | Unit | Description |
|----------------------------------------------------------|-------------|------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
| `RunID` | integer | - | Run identifier (ensure these are all unique values for every row) |
| `ModelName` | string | - | Name for this analysis. Name is used to create output directory with results |
| `NodeFile` | string | - | File name with node input data, including .csv extension |
| `SegmentFile` | string | - | File name with segment input data, including .csv extension |
| `phi` | double | rad | Soil peak angle of internal friction (used to calculate lateral resistance) |
| `phi_cs` | double | rad | Soil critical state angle (used to calculate effect of soil confinement on additional shear resistance) |
| `delta` | double | rad | Interface friction angle between soil and root |
| `gamma` | double | N/mm3 | Soil effective unit weight |
| `K` | double | - | Coefficient of lateral earth pressure |
| `sigmav0` | double | MPa | Vertical effective stress at the soil surface |
| `shearplane_param` | double | mm | Shear plane thickness parameter (b_sh) in paper |
| `shearplane_depth` | double | mm | Shear plane depth, at Y=0 |
| `shearplane_angle` | double | rad | Shear plane angle |
| `n_node` | integer | - | Number of nodes per segment |
| `n_step_max` | integer | - | Maximum number of displacement steps |
| `n_iter_max` | integer | - | maximum number of iterations (to find solution for a single displacement step) |
| `uext_inc0` | double | mm | Initial shear displacement step size |
| `uext_max` | double | mm | Maximum shear displacement (analysis stops once this value is reached) |
| `uext_factor_increase` | double | - | Multiplaction factor to automatically increase displacement step size if suitable |
| `uext_factor_decrease` | double | - | Multiplaction factor to decrease displacement step if solution if solution not found |
| `uext_incmin` | double | mm | Minimum displacement step size (if lower required, analysis stops) |
| `uext_incmax` | double | mm | Maximum displacement step size |
| `solve_tolerance` | double | - | Solver tolerance of &lt;scipy.integrate.solve_bvp&gt; solver in Python |
| `SaveStep` | integer | - | Save full root and soil data every &lt;SaveStep&gt; displacement step |`

## Root geometry
The root architecture is described by two seperate files:
* Nodal data   - This defined the position and displacement boundary conditions of individual nodes (points) in the geometry
* Segment data - This defines the connections between two nodes. For each connection, segment properties such as stiffness of diameter are given
Theoretically, there is no limit to the amount of segments, although it can be expected that the code will become less efficient with many segments (as the number of differential equations to solve increases)

### Node input data format
Required columns:

| Parameter | Type | Unit | Description |
|----------------------------------------------------------|-------------|------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------|
| `NodeID` | integer | - | Node identifier (ensure these are all unique values) |
| `X` | double | mm | Global X-position of node |
| `Y` | double | mm | Global Y-position of node |
| `bound_X` | integer | - | X-displacement constrained (value=1) or free (0) |
| `bound_Y` | integer | - | Y-displacement constrained (value=1) or free (0) |
| `bound_Theta` | integer | - | rotation constrained (value=1) or free (0) |

### Segment input data format
Required columns:

| Parameter | Type | Unit | Description |
|----------------------------------------------------------|-------------|------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------|
| `SegmentID` | integer | - | Segment identifier (ensure these are all unique values) |
| `NodeID1` | double | - | Node identifier of node at start of segment |
| `NodeID2` | double | - | Node identifier of node at end of segment |
| `d` | double | mm | Root diameter |
| `Et1` | double | MPa | Coefficient in root tensile stiffness (\( \xi_1 \) in paper) |
| `Et2` | double | - | Coefficient in root tensile stiffness (\xi_2 in paper) |
| `Et3` | double | MPa | Coefficient in root tensile stiffness (\xi_3 in paper) |
| `Eb1` | double | MPa | Coefficient in root bending stiffness (\xi_1 in paper) |
| `Eb2` | double | - | Coefficient in root bending stiffness (\xi_2 in paper) |
| `Eb3` | double | MPa | Coefficient in root bending stiffness (\xi_3 in paper) |

# General remarks on how the code works
* The program loads all input data in input parameter file specified
Every line in this file represents a single analysis (i.e. a single direct shear tests)
Required parameters are given above, as are node and segment parameters for the root architecture
Subsequently, the code runs through every run in sequence
* The set of differential equations and boundary conditions is solved for every displacement step. Initial step size is `uext_inc0`
The initial guess required is based on the last solution
When no solutions can be found, the displacement step is decreased (multiplying current step times `uext_factor_decrease`) to increase the chance that the system will be solved.
* If the displacement step becomes too small (smaller than `uext_incmin`), the analysis for this run will automatically stop.
The displacement step is multiplied by `uext_factor_increase` when steady solutions are found systematically, but will never be larger than `uext_incmax`.
The analysis stops when no solution can be found within acceptable step sizes, or when the final displacement `uext_max` is reached
* The set of differential equations is written in a set of first order differential equations, as is a requirement for the `solve_bvp` solver in SciPy. All differential equations are normalised by segment length.
* The order of variables is: theta, dtheta/ds. d^2theta/ds^2, eps, u, w. When there are multiple segments, these are added in sequence of occurance in the `dsegm` segment data dictionary (so first all variables for the the first segments, than all those for the second etc.)

# Outputs
Depending on the output settings in `Main.py`, output will be generated and saved. 
When `data_save` equals `True`, the code prints output `.csv` files to a new directory `Results`. 
Subfolders will be created based on the `ModelName` and `RunID` given in the input parameter file.
Data will be saved for every `SaveStep` displacement step for both soil and root, as well as summary data for the whole analyses.

## Output - Root
For every displacement step, a `.csv` is created containing the following fields:

| Parameter | Type | Unit | Description |
|----------------------------------------------------------|-------------|------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------|
| `s` | double | mm | Local coordinate along the (displaced) root axis |
| `t` | double | rad | Root rotation (root coordinate system) |
| `dt` | double | rad/mm | Root curvature (dt/ds) |
| `ddt` | double | rad/mm2 | Derivative of root curvature (d^2t/ds^2) |
| `eps` | double | - | Root axial strain |
| `u` | double | mm | Root axial displacement (root coordinate system) |
| `w` | double | mm | Root lateral displacement (root coordinate system) |
| `x` | double | mm | Root x-position (root coordinate system) |
| `y` | double | mm | Root y-position (root coordinate system) |
| `Xr` | double | mm | Root X-position (global coordinate system) |
| `Yr` | double | mm | Root Y-position (global coordinate system) |
| `umob` | double | mm | Relative axial soil-root displacement (delta u in paper) |
| `wmob` | double | mm | Relative lateral soil-root displacement (delta w in paper) |
| `SegmentID` | integer | - | Segment identifier |

## Output - Soil
For every displacement step, a `.csv` is created containing the following fields:

| Parameter | Type | Unit | Description |
|----------------------------------------------------------|-------------|------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------|
| `Xs` | double | mm | Soil X-position (global coordinate system) |
| `Ys` | double | mm | Soil Y-position (global coordinate system) |
| `SegmentID` | integer | - | Segment identifier |

## Output - Summary
Three summary data files are outputted
* `ModelName` + `_Results_ShearForceDisp.csv`:

| Parameter | Type | Unit | Description |
|----------------------------------------------------------|-------------|------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------|
| `StepID` | integer | mm | Displacement step identifier |
| `uext` | double | mm | Shear displacement |
| `F_total` | double | N | Total root-reinforcement in middle of shear plane |

* `ModelName` + `_Results_ShearReinforcement.csv`:

| Parameter | Type | Unit | Description |
|----------------------------------------------------------|-------------|------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------|
| `StepID` | integer | mm | Displacement step identifier |
| `SegmentID` | integer | - | Segment identifier |
| `F_direct` | double | N | Direct contribution of force in root to shear resistance (force component along shear displacement direction) |
| `F_comfine` | double | N | Additional normal force on shear plane due to roots (positive force compresses soil further) |
| `F_total` | double | N | Total root-reinforcement in middle of shear plane |

* `ModelName` + `_Results_Stresses.csv`:

| Parameter | Type | Unit | Description |
|----------------------------------------------------------|-------------|------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------|
| `StepID` | integer | mm | Displacement step identifier |
| `SegmentID` | integer | - | Segment identifier |
| `sig_ax_shearplane` | double | MPa | Axial stress at shear plane |
| `sig_be_shearplane` | double | MPa | Bending stress at shear plane |
| `sig_sh_shearplane` | double | MPa | Max. shear stress at shear plane |
| `sig_ax_max` | double | MPa | Max axial stress in segment |
| `sig_be_max` | double | MPa | Max bending stress in segment |
| `sig_sh_max` | double | MPa | Max shear stress in segment |
| `sig_fi_max` | double | MPa | Max stress in ultimate fibre in segment |
| `sig_fi_max_ax` | double | MPa | Axial stress component of `sig_fi_max` |
| `sig_fi_max_be` | double | MPa | Bending stress component of `sig_fi_max` |
| `X_sig_ax_max` | double | mm | Global (displaced) X-position where `sig_ax_max` occurs |
| `X_sig_be_max` | double | mm | Global (displaced) X-position where `sig_be_max` occurs |
| `X_sig_sh_max` | double | mm | Global (displaced) X-position where `sig_sh_max` occurs |
| `X_sig_fi_max` | double | mm | Global (displaced) X-position where `sig_fi_max` occurs |


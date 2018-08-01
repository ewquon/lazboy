from __future__ import print_function
import os
import yaml

mypath = os.path.dirname(os.path.realpath(__file__))

def get_templates():
    dpath = os.path.join(mypath, 'simulation_templates')
    listing = os.listdir(dpath)
    templates = { os.path.splitext(name)[0]: os.path.join(dpath,name)
                    for name in listing
                    if os.path.isfile(os.path.join(dpath,name))
                        and name.endswith('.yaml') }

    listing = os.listdir(os.getcwd())
    print(listing)
    custom_configs = { os.path.splitext(name)[0]: os.path.join(os.getcwd(),name)
                        for name in listing
                        if os.path.isfile(name) and name.endswith('.yaml') }

    return templates, custom_configs

def read_template(fpath):
    #fpath = os.path.join(mypath, 'simulation_templates', name+'.yaml')
    print('Loaded template: '+fpath)
    with open(fpath,'r') as f:
        params = yaml.load(f)
    return params

def save_template(name):
    fpath = os.path.join(mypath, 'simulation_templates', name)
    if os.path.isfile(fpath):
        print('Template '+name+' already exists! No file was written.')
        return

yaml_template = """#
# SOWFA precursor template generated by sowfa_precursor_setup.py
# github.com/ewquon/lazboy
#

# domain controls
xMin: {xMin:f}  # Minimum x-extent of domain (m).
yMin: {yMin:f}  # Minimum y-extent of domain (m).
zMin: {zMin:f}  # Minimum z-extent of domain (m).
xMax: {xMax:f}  # Maximum x-extent of domain (m).
yMax: {yMax:f}  # Maximum y-extent of domain (m).
zMax: {zMax:f}  # Maximum z-extent of domain (m).
nx: {nx:d}  # Number of cells in x-direction.
ny: {ny:d}  # Number of cells in y-direction.
nz: {nz:d}  # Number of cells in z-direction.


# decomposition controls
PPN: {PPN:d}  # Processors per node
nCores: {nCores:d}  # Number of cores on which to run this case.
decompType: {decompType:s}  # Decomposition algorithm.  "simple" and "scotch" are good choices.
decompOrder: [{decompOrder[0]:d}, {decompOrder[1]:d}, {decompOrder[2]:d}] # Order of the decomposition number of partitions in (x y z)-directions.


# general conditions
TRef: {TRef:f}  # Reference potential temperature (K).
coriolis: {coriolis}
latitude: {latitude:f}  # Latitude on the Earth of the site (deg).
EarthPeriod: {EarthPeriod:f}  # Earth's rotation period (hr).


# atmosphere controls
U0Mag: {U0Mag:f}  # Initial condition for wind speed (m/s).
dir: {dir:f}  # Initial condition for wind direction (deg).
windHeight: {windHeight:f}  # Height at which to drive mean wind to U0Mag/dir (m).
p_rgh0: {p_rgh0:f}  # Initial pressure (minus the hydrostatic variation and normalized by density) (m^2/s^2).
nuSgs0: {nuSgs0:f}  # Initial SGS viscosity (m^2/s).
k0: {k0:f}  # Initial SGS turbulent kinetic energy (m^2/s^2).
kappat0: {kappat0:f}  # Initial SGS temperature diffusivity (m^2/s).
TGradUpper: {TGradUpper:f}  # Potential temperature gradient above the strong inversion (K/m).
zInversion: {zInversion:f}  # Height of the middle of the initial strong capping inversion (m).
inversionWidth: {inversionWidth:g}  # Vertical width of the intial strong capping inversion (m).
TBottom: {TBottom:f}  # Initial potential temperature at bottom of strong capping inversion (K).
TTop: {TTop:f}  # Initial potential temperature at top of strong capping inversion (K).

sourceType: {sourceType:s}
idealProfile: {idealProfile}
alpha: {alpha:f}  # Shear exponent.
veer: {veer:f}  # Veer (deg).
profileTable: {profileTable}


# surface controls
surfaceBCType: {surfaceBCType:s}
qwall: {qwall:f}  # Temperature flux at wall (modify the z-value).  A negative value is flux into domain (K-m/s).
Rwall: {Rwall}  # Initial wall shear stress (m^2/s^2).
kappa: {kappa:f}  # von Karman constant.
z0: {z0:g}  # Surface roughness (m).
heatingRate: {heatingRate:g}  # Surface temperature change rate (when not directly setting temperature flux) (K/s).


# advanced settings

# -surface conditions
wallModelAverageType: {wallModelAverageType:s}  # Treat surface stress wall model locally ("local") or with planar averaging ("planarAverage").
betaM: {betaM:f}  # Monin-Obukhov wall shear stress model constant.
gammaM: {gammaM:f}  # Monin-Obukhov wall shear stress model constant.
betaH: {betaH:f}  # Monin-Obukhov wall temperature flux model constant.
gammaH: {gammaH:f}  # Monin-Obukhov wall temperature flux model constant.
alphaH: {alphaH:f}  # Monin-Obukhov wall temperature flux model constant.

# -planar averaging and source term statistics options.
statisticsOn: {statisticsOn}  # Gather planar-averaged flow statistics.
statisticsFrequency: {statisticsFrequency:d}  # Frequency in time steps of statistics gathering.

# -transport properties
Pr: {Pr:f}  # Molecular Prandtl number.
Prt: {Prt:f}  # Turbulent Prandtl number.
nu: {nu:g}  # Molecular viscosity (m^2/s).

# -SGS model inputs
LESModel: {LESModel:s}  # SGS model selection.
ce: {ce:g}  # SGS model constant.
ck: {ck:g}  # SGS model constant.
"""


setUp_template = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.x                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/

// Domain size and number of cells.
xMin                 {xMin:f};                         // Minimum x-extent of domain (m).
yMin                 {yMin:f};                         // Minimum y-extent of domain (m).
zMin                 {zMin:f};                         // Minimum z-extent of domain (m).
xMax                 {xMax:f};                      // Maximum x-extent of domain (m).
yMax                 {yMax:f};                      // Maximum y-extent of domain (m).
zMax                 {zMax:f};                      // Maximum z-extent of domain (m).
nx                   {nx:d};                         // Number of cells in x-direction.
ny                   {ny:d};                         // Number of cells in y-direction.
nz                   {nz:d};                         // Number of cells in z-direction.




// Number of cores and domain decomposition information.
nCores               {nCores:d};                         // Number of cores on which to run this case.
decompType           {decompType:s};                      // Decomposition algorithm.  "simple" and "scotch" are good choices.
decompOrder          ({decompOrder[0]:d} {decompOrder[1]:d} {decompOrder[2]:d});                   // Order of the decomposition number of partitions in (x y z)-directions.




// Planar averaging and source term statistics options.
statisticsOn         {statisticsOn};                        // Gather planar-averaged flow statistics.
statisticsFrequency  {statisticsFrequency:d};                           // Frequency in time steps of statistics gathering.




// Initial values for the variables.
// Note that U and T get overwritten if setFieldsABL is called.
U0Mag                {U0Mag:f};                         // Initial condition for wind speed (m/s).
dir                  {dir:f};                       // Initial condition for wind direction (deg).
windHeight           {windHeight:f};                       // Height at which to drive mean wind to U0Mag/dir (m).
p_rgh0               {p_rgh0:f};                         // Initial pressure (minus the hydrostatic variation and normalized by density) (m^2/s^2).
nuSgs0               {nuSgs0:f};                         // Initial SGS viscosity (m^2/s).
k0                   {k0:f};                         // Initial SGS turbulent kinetic energy (m^2/s^2).
kappat0              {kappat0:f};                         // Initial SGS temperature diffusivity (m^2/s).
TGradUpper           {TGradUpper:f};                       // Potential temperature gradient above the strong inversion (K/m).
zInversion           {zInversion:f};                       // Height of the middle of the initial strong capping inversion (m).
inversionWidth       {inversionWidth:f};                       // Vertical width of the intial strong capping inversion (m).
TBottom              {TBottom:f};                       // Initial potential temperature at bottom of strong capping inversion (K).
TTop                 {TTop:f};                       // Initial potential temperature at top of strong capping inversion (K).




// General conditions and parameters.
Pr                   {Pr:f};                         // Molecular Prandtl number.
Prt                  {Prt:f};                  // Turbulent Prandtl number.
nu                   {nu:g};                      // Molecular viscosity (m^2/s).
TRef                 {TRef:f};                       // Reference potential temperature (K).
latitude             {latitude:f};                        // Latitude on the Earth of the site (deg).
EarthPeriod          {EarthPeriod:f};                        // Earth's rotation period (hr).




// SGS model inputs.
LESModel             {LESModel:s};                // SGS model selection.
ce                   {ce:f};                        // SGS model constant.
ck                   {ck:f};                      // SGS model constant.




// Surface conditions.
qwall               (0.0 0.0 {qwall:f});                // Temperature flux at wall (modify the z-value).  A negative value is flux into domain (K-m/s).
Rwall               ({Rwall[0]:f} {Rwall[1]:f} {Rwall[2]:f} {Rwall[3]:f} {Rwall[4]:f} {Rwall[5]:f});    // Initial wall shear stress (m^2/s^2).
kappa                {kappa:f};                         // von Karman constant.
z0                   {z0:f};                        // Surface roughness (m).
wallModelAverageType "{wallModelAverageType:s}";             // Treat surface stress wall model locally ("local") or with planar averaging ("planarAverage").
betaM                {betaM:f};                        // Monin-Obukhov wall shear stress model constant.
gammaM               {gammaM:f};                         // Monin-Obukhov wall shear stress model constant.
betaH                {betaH:f};                         // Monin-Obukhov wall temperature flux model constant.
gammaH               {gammaH:f};                         // Monin-Obukhov wall temperature flux model constant.
alphaH               {alphaH:f};                         // Monin-Obukhov wall temperature flux model constant.
heatingRate          {heatingRate:g};                         // Surface temperature change rate (when not directly setting temperature flux) (K/s).




#inputMode           merge

// ************************************************************************* //
"""

runscript_preprocess_template = """#!/bin/bash
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn={PPN:d}
#PBS -l feature={PPN:d}core
#PBS -A {allocation:s}
#PBS -q batch
#PBS -m abe
#PBS -M {email:s}


# User Input.
OpenFOAMversion=2.4.x           # OpenFOAM version
startTime=0                     # Start time
updateBCType=0                  # Boolean for whether or not the boundary condition types will be updated over 
                                #    what is in the initial conditions files. Leave it 0 for precursors
inflowDir='cyclic'              # For inflow/outflow cases, specify the inflow direction.  Choices are 'west',
                                #   'east', 'south', 'west', 'southWest', 'northWest', 'southEast', and
                                #   'northEast'.  There is a 'cyclic' option too in case you need to change
                                #    back to cyclic or simply update the number of boundary face entries.
parallel=0                      # Boolean for whether or not the preprocessing is run in parallel.
cores=1                         # Enter the number of cores you will preprocess on.

refinementLevels=0              # If you want to refine the mesh locally for any reason, or if you are making
                                # a uniform resolution mesh that is so large that you need to build it in serial
                                # at lower resolution and then globally refine, set the number of refinement
                                # levels here.  See the refineMeshSerial and refineMeshParallel functions to 
                                # understand what they do.  The serial version runs topoSet and refineHexMesh, 
                                # so you need to provide system/topoSetDict.local.N files where N corresponds
                                # to the refinement level (i.e., if you are doing nested local refinement boxes.
                                # In most cases, though, you probably will not be refining, so keep this set to 
                                # 0.


# Define some functions for mesh refinement.
# Local refinement performed on one core.
refineMeshLocal()
{{
   i=$1
   while [ $i -ge 1 ]
   do
      echo "   -Performing level $i local refinement with topoSet/refineHexMesh"
      echo "      *selecting cells to refine..."
      topoSet -dict system/topoSetDict.local.$i > log.topoSet.local.$i 2>&1

      echo "      *refining cells..."
      refineHexMesh local -overwrite > log.refineHexMesh.local.$i 2>&1

      let i=i-1
   done
}}

# Global refinement performed in parallel.
refineMeshGlobal()
{{
   i=1
   while [ $i -le $1 ]
   do
      echo "   -Performing level $i global refinement with refineMesh"
      echo "      *refining cells..."
      mpirun -np $cores refineMesh -all -parallel -overwrite > log.refineMesh.global.$i 2>&1

      let i=i+1
   done
}}


# If running in parallel, cd to job launch directory
if [ $parallel -eq 1 ]
   then
   cd $PBS_O_WORKDIR
fi


# Source the bash profile and then call the appropriate OpenFOAM version function
# so that all the modules and environment variables get set.
source $HOME/.bash_profile
OpenFOAM-$OpenFOAMversion


# Copy the controlDict.1 (assuming this is the one the actual solver will start
# out with) to controlDict.  OpenFOAM reads "controlDict", not the numbered ones.
cp system/controlDict.1 system/controlDict


# Copy the "clean" .original initial fields to a working copy.  OpenFOAM does not
# read the ".original" initial fields--that's why they remain clean.
rm -rf $startTime
cp -rf $startTime.original $startTime


# Build the mesh.
cp constant/polyMesh/blockMeshDict ./
rm -rf constant/polyMesh/*
mv ./blockMeshDict constant/polyMesh
blockMesh > log.blockMesh 2>&1


# The initial fields come from the precursor which is periodic on all sides.  The turbine
# case has inflow and outflow.  Call the changeDictionary utility to make the south and
# north sides inflow and outflow.
if [ $updateBCType -eq 1 ]
   then
   changeDictionary -dict system/changeDictionaryDict.updateBCs.$inflowDir -time $startTime -enableFunctionEntries > log.changeDictionary.updateBCs.$inflowDir.1 2>&1
fi


# Do serial local refinement
refineMeshLocal $refinementLevels


# If running in parallel from this point forward, then do the following:
if [ $parallel -eq 1 ]
   then
   # Decompose the mesh and solution files (serial)
   decomposePar -cellDist -force > log.decomposePar 2>&1

   # Check the mesh
   mpirun -np $cores checkMesh -parallel > log.checkMesh.1 2>&1

   # Perform global refinement to desired resolution.
   refineMeshGlobal $refinementLevels

   # The mesh got globally refined, but the solution file did not, so
   # the boundary fields may not have the correct number of entries.
   # Use the changeDictionary utility to overwrite the spatially varying
   # boundary data to a uniform single value.
   if [ $updateBCType -eq 1 ]
      then
      mpirun -np $cores changeDictionary -dict system/changeDictionaryDict.updateBCs.$inflowDir -time $startTime -enableFunctionEntries -parallel > log.changeDictionary.updateBCs.$inflowDir.1 2>&1
   fi

   # Renumber the mesh for better matrix solver performance.
   mpirun -np $cores renumberMesh -parallel -overwrite > log.renumberMesh 2>&1

   # Do one last check on the mesh.
   mpirun -np $cores checkMesh -parallel > log.checkMesh.3 2>&1


# Otherwise, run in serial as follows:
else
   # Renumber the mesh.
   echo "   -Renumbering the mesh with renumberMesh..."
   renumberMesh -overwrite > log.renumberMesh 2>&1

   # Decompose the mesh and solution files (serial)
   echo "   -Decomposing the domain with decomposePar..."
   decomposePar -cellDist -force > log.decomposePar 2>&1

   # Check the mesh.
   echo "   -Checking the mesh with checkMesh..."
   checkMesh > log.checkMesh.1 2>&1
fi
"""

runscript_solve_template = """#!/bin/bash
#PBS -N {casename:s}
#PBS -l walltime=48:00:00
#PBS -l nodes={nNodes:d}:ppn={PPN:d}
#PBS -l feature={PPN:d}core
#PBS -A {allocation:s}
#PBS -q batch
#PBS -m abe
#PBS -M {email:s}

source $HOME/.bash_profile
OpenFOAM-2.4.x
module list

cd $PBS_O_WORKDIR

cores={nCores:d}

initializer=setFieldsABL
solver=ABLSolver
runNumber=1
startTime=0

cp system/controlDict.$runNumber system/controlDict

echo "Starting OpenFOAM job at: " $(date)
echo "using " $cores " cores"

# Run the flow field initializer (parallel)
if [ $runNumber -eq 1 ] 
   then
   mpirun -np $cores $initializer -parallel > log.$runNumber.$initializer 2>&1
fi

# Run the solver (parallel)
mpirun -np $cores $solver -parallel > log.$runNumber.$solver 2>&1

echo "Ending OpenFOAM job at: " $(date)
"""

setFieldsABLDict_template = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.0.x                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/

FoamFile
{{
    version         2.0;
    format          ascii;

    root            "";
    case            "";
    instance        "";
    local           "";

    class           dictionary;
    object          setFieldsABLDict;
}}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include        "../setUp"


// Extents of the domain.
xMin                      $xMin;
yMin                      $yMin;
zMin                      $zMin;

xMax                      $xMax;
yMax                      $yMax;
zMax                      $zMax;

zRef                      $zMax;

// Specify if distance from wall should be used as z.
useWallDistZ              false;
scaleVelocityWithHeight   false;


// Specify how to initialze the base velocity and temperature profile.
velocityInitType          "{velocityInitType:s}";
temperatureInitType       "{temperatureInitType:s}";

// Maximum perturbation of streamwise/spanwise flow near surface.
deltaU                    0.25;
deltaV                    0.25;

// Total periods of perturbations in streamwise/spanwise in the domain.
Uperiods                  12.0;
Vperiods                  12.0;

// Percentage of domain height (zMax) where peak in perturbation 
// magnitude occurs.
zPeak                     0.015;

// Initial height of the center of the capping inversion.
zInversion                $zInversion;

// Width of the capping inversion.
widthInversion            $inversionWidth;

// Potential temperature at the bottom and top of the capping inversion.
Tbottom                   $TBottom;
Ttop                      $TTop;

// Maximum temperature fluctuation size below capping inversion.
TPrimeScale               0.0;

// Height rate of change of potential temperature above the inversion.
dTdz                      $TGradUpper;

// Geostrophic wind speed magnitude.
Ug                        $U0Mag;

// Geostrophic wind speed direction.
UgDir                     $dir;

// Aerodynamic roughness height of surface.
z0                        $z0;

// von Karman constant.
kappa                     $kappa;

// Vertical profile table.
profileTable
(
{profileTable:s}
);

// Update internal field.
updateInternalFields      true;

// Update boundary field.
updateBoundaryFields      false;

// ************************************************************************* //
"""

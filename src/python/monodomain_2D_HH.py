#!/usr/bin/env python

#DOC-START imports
import sys, os, math

# setting random seed to recreate the results for parallel runs
import random
random.seed(100)

# Intialise OpenCMISS
from opencmiss.iron import iron
#DOC-END imports

# Set problem parameters
#DOC-START parameters
# 2D domain size
height = 0.05 #cm
width = 0.1 #cm
numberOfXElements = 25
numberOfYElements = 13

# Materials parameters
Am = 1936 #cm^-1
Cm = 1.4 #uF.cm^2
conductivity = 1.0 #mS.cm^-1
# Simulation parameters
# stimValue is directly passed on to the cellml equations
stimValue = -5000.0 #uA.cm^-2
stimStop = 0.2 #ms
timeStop = 0.7 #ms
odeTimeStep = 0.00001 #ms
pdeTimeStep = 0.001 #ms
outputFrequency = 10
#DOC-END parameters

#Setup field number handles
coordinateSystemUserNumber = 1
regionUserNumber = 1
basisUserNumber = 1
pressureBasisUserNumber = 2
generatedMeshUserNumber = 1
meshUserNumber = 1
cellMLUserNumber = 1
decompositionUserNumber = 1
equationsSetUserNumber = 1
problemUserNumber = 1
#Mesh component numbers
linearMeshComponentNumber = 1
#Fields
geometricFieldUserNumber = 1
fibreFieldUserNumber = 2
dependentFieldUserNumber = 3
materialsFieldUserNumber = 4
equationsSetFieldUserNumber = 5
cellMLModelsFieldUserNumber = 6
cellMLStateFieldUserNumber = 7
cellMLParametersFieldUserNumber = 8
cellMLIntermediateFieldUserNumber = 9

#DOC-START parallel information
# Get the number of computational nodes and this computational node number
computationEnvironment = iron.ComputationEnvironment()
numberOfComputationalNodes = computationEnvironment.NumberOfWorldNodesGet()
computationalNodeNumber = computationEnvironment.WorldNodeNumberGet()
#DOC-END parallel information

#DOC-START initialisation
# Create a 2D rectangular cartesian coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.DimensionSet(2)
coordinateSystem.CreateFinish()

# Create a region and assign the coordinate system to the region
region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.LabelSet("Region")
region.coordinateSystem = coordinateSystem
region.CreateFinish()
#DOC-END initialisation

#DOC-START basis
# Define a bilinear Lagrange basis
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
basis.numberOfXi = 2
basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*2
basis.quadratureNumberOfGaussXi = [3]*2
basis.CreateFinish()
#DOC-END basis

#DOC-START generated mesh
# Create a generated mesh
generatedMesh = iron.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber,region)
generatedMesh.type = iron.GeneratedMeshTypes.REGULAR
generatedMesh.basis = [basis]
generatedMesh.extent = [width,height]
generatedMesh.numberOfElements = [numberOfXElements,numberOfYElements]

mesh = iron.Mesh()
generatedMesh.CreateFinish(meshUserNumber,mesh)
#DOC-END generated mesh

#DOC-START decomposition
# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()
#DOC-END decomposition

#DOC-START geometry
# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber, region)
geometricField.meshDecomposition = decomposition
geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(iron.FieldVariableTypes.U, "coordinates")
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, linearMeshComponentNumber)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, linearMeshComponentNumber)
geometricField.CreateFinish()

# Set geometry from the generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)
#DOC-END geometry

#DOC-START equations set
# Create the equations_set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.BIOELECTRICS,
        iron.EquationsSetTypes.MONODOMAIN_EQUATION,
        iron.EquationsSetSubtypes.NONE]
equationsSet.CreateStart(equationsSetUserNumber, region, geometricField,
        equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()
#DOC-END equations set

#DOC-START equations set fields
# Create the dependent Field
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
equationsSet.DependentCreateFinish()

# Create the materials Field
materialsField = iron.Field()
equationsSet.MaterialsCreateStart(materialsFieldUserNumber, materialsField)
equationsSet.MaterialsCreateFinish()

# Set the materials values
# Set Am
materialsField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,Am)
# Set Cm
materialsField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,Cm)
# Set conductivity
materialsField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,conductivity)
materialsField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,4,conductivity)
#DOC-END equations set fields

# Read the cellml file either as an argument (useful for testing) or hardcoded text.
if len(sys.argv) > 1:
	cellmlFile = sys.argv[1]
else:
	cellmlFile = './HodgkinHuxley1952.cellml'

#DOC-START create cellml environment
# Create the CellML environment
cellML = iron.CellML()
cellML.CreateStart(cellMLUserNumber, region)
# Import the cell model from a file
cellModel = cellML.ModelImport(cellmlFile)
#DOC-END create cellml environment

#DOC-START flag variables
# Now we have imported the model we are able to specify which variables from the model we want to set from openCMISS
cellML.VariableSetAsKnown(cellModel, "membrane/i_Stim")
cellML.VariableSetAsKnown(cellModel, "membrane/Cm")
# and variables to get from the CellML
cellML.VariableSetAsWanted(cellModel, "membrane/i_K")
#DOC-END flag variables

#DOC-START create cellml finish
cellML.CreateFinish()
#DOC-END create cellml finish

#DOC-START map Vm components
# Start the creation of CellML <--> OpenCMISS field maps
cellML.FieldMapsCreateStart()
#Now we can set up the field variable component <--> CellML model variable mappings.
#Map Vm
cellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U,1, iron.FieldParameterSetTypes.VALUES,cellModel,"membrane/V", iron.FieldParameterSetTypes.VALUES)
cellML.CreateCellMLToFieldMap(cellModel,"membrane/V", iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#Finish the creation of CellML <--> OpenCMISS field maps
cellML.FieldMapsCreateFinish()

# Set the initial Vm values
dependentField.ComponentValuesInitialise(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,-75)

#DOC-END map Vm components

#DOC-START define CellML models field
#Create the CellML models field
cellMLModelsField = iron.Field()
cellML.ModelsFieldCreateStart(cellMLModelsFieldUserNumber, cellMLModelsField)
cellML.ModelsFieldCreateFinish()
#DOC-END define CellML models field

#DOC-START define CellML state field
#Create the CellML state field
cellMLStateField = iron.Field()
cellML.StateFieldCreateStart(cellMLStateFieldUserNumber, cellMLStateField)
cellML.StateFieldCreateFinish()
#DOC-END define CellML state field

#DOC-START define CellML parameters and intermediate fields
#Create the CellML parameters field
cellMLParametersField = iron.Field()
cellML.ParametersFieldCreateStart(cellMLParametersFieldUserNumber, cellMLParametersField)
cellML.ParametersFieldCreateFinish()

# Set the initial Cm values
cmComponent = cellML.FieldComponentGet(cellModel, iron.CellMLFieldTypes.PARAMETERS, "membrane/Cm")
cellMLParametersField.ComponentValuesInitialise(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, cmComponent,Cm)

#  Create the CellML intermediate field
cellMLIntermediateField = iron.Field()
cellML.IntermediateFieldCreateStart(cellMLIntermediateFieldUserNumber, cellMLIntermediateField)
cellML.IntermediateFieldCreateFinish()
#DOC-END define CellML parameters and intermediate fields

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

# Find the domains of the first and last nodes
firstNodeNumber = 1
lastNodeNumber = (numberOfXElements+1)*(numberOfYElements+1)
firstNodeDomain = decomposition.NodeDomainGet(firstNodeNumber, 1)
lastNodeDomain = decomposition.NodeDomainGet(lastNodeNumber, 1)

# Set the stimulus on half the bottom nodes
stimComponent = cellML.FieldComponentGet(cellModel, iron.CellMLFieldTypes.PARAMETERS, "membrane/i_Stim")
for node in range(1,(numberOfXElements + 1)/2 + 1):
    nodeDomain = decomposition.NodeDomainGet(node,1)
    if nodeDomain == computationalNodeNumber:
        cellMLParametersField.ParameterSetUpdateNode(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, node, stimComponent, stimValue)

#DOC-START define monodomain problem
#Define the problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.BIOELECTRICS,
    iron.ProblemTypes.MONODOMAIN_EQUATION,
    iron.ProblemSubtypes.MONODOMAIN_GUDUNOV_SPLIT]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()
#DOC-END define monodomain problem

#Create the problem control loop
problem.ControlLoopCreateStart()
controlLoop = iron.ControlLoop()
problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],controlLoop)
controlLoop.TimesSet(0.0,stimStop,pdeTimeStep)
controlLoop.OutputTypeSet(iron.ControlLoopOutputTypes.TIMING)
controlLoop.TimeOutputSet(outputFrequency)
problem.ControlLoopCreateFinish()

#Create the problem solvers
daeSolver = iron.Solver()
dynamicSolver = iron.Solver()
problem.SolversCreateStart()
# Get the first DAE solver
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,daeSolver)
daeSolver.DAETimeStepSet(odeTimeStep)
daeSolver.OutputTypeSet(iron.SolverOutputTypes.NONE)
# Get the second dynamic solver for the parabolic problem
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,dynamicSolver)
dynamicSolver.OutputTypeSet(iron.SolverOutputTypes.NONE)
problem.SolversCreateFinish()

#DOC-START define CellML solver
#Create the problem solver CellML equations
cellMLEquations = iron.CellMLEquations()
problem.CellMLEquationsCreateStart()
daeSolver.CellMLEquationsGet(cellMLEquations)
cellmlIndex = cellMLEquations.CellMLAdd(cellML)
problem.CellMLEquationsCreateFinish()
#DOC-END define CellML solver

#Create the problem solver PDE equations
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
dynamicSolver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Prescribe any boundary conditions
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
solverEquations.BoundaryConditionsCreateFinish()

# Solve the problem until stimStop
problem.Solve()

# Now turn the stimulus off
for node in range(1,(numberOfXElements + 1)/2 + 1):
#for node in range(1,numberOfXElements+2):
    nodeDomain = decomposition.NodeDomainGet(node,1)
    if nodeDomain == computationalNodeNumber:
        cellMLParametersField.ParameterSetUpdateNode(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, node, stimComponent, 0.0)

#Set the time loop from stimStop to timeStop
controlLoop.TimesSet(stimStop,timeStop,pdeTimeStep)

# Now solve the problem from stim stop until time stop
problem.Solve()

# Export the results, here we export them as standard exnode, exelem files
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("Monodomain","FORTRAN")
fields.ElementsExport("Monodomain","FORTRAN")
fields.Finalise()

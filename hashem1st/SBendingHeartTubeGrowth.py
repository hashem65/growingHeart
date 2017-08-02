#!/usr/bin/env python

#> \file
#> \author Chris Bradley
#> \brief This is an example script to illustrate embryonic S-bending in a heart tube using OpenCMISS calls.
#>
#> \section LICENSE
#>
#> Version: MPL 1.1/GPL 2.0/LGPL 2.1
#>
#> The contents of this file are subject to the Mozilla Public License
#> Version 1.1 (the "License"); you may not use this file except in
#> compliance with the License. You may obtain a copy of the License at
#> http://www.mozilla.org/MPL/
#>
#> Software distributed under the License is distributed on an "AS IS"
#> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#> License for the specific language governing rights and limitations
#> under the License.
#>
#> The Original Code is OpenCMISS
#>
#> The Initial Developer of the Original Code is University of Auckland,
#> Auckland, New Zealand and University of Oxford, Oxford, United
#> Kingdom. Portions created by the University of Auckland and University
#> of Oxford are Copyright (C) 2007 by the University of Auckland and
#> the University of Oxford. All Rights Reserved.
#>
#> Contributor(s): 
#>
#> Alternatively, the contents of this file may be used under the terms of
#> either the GNU General Public License Version 2 or later (the "GPL"), or
#> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
#> in which case the provisions of the GPL or the LGPL are applicable instead
#> of those above. if you wish to allow use of your version of this file only
#> under the terms of either the GPL or the LGPL, and not to allow others to
#> use your version of this file under the terms of the MPL, indicate your
#> decision by deleting the provisions above and replace them with the notice
#> and other provisions required by the GPL or the LGPL. if you do not delete
#> the provisions above, a recipient may use your version of this file under
#> the terms of any one of the MPL, the GPL or the LGPL.
#>

#> Main script
# Add Python bindings directory to PATH
import sys, os
import math
import numpy as np

# Intialise OpenCMISS
from opencmiss.iron import iron

#iron.OutputSetOn("Testing")

# Set problem parameters
growthModel = 1 # Type of growth. 1 - volumetric; 2 - stress based; 3 - strain based.
isotropic = False # True if the problem is isotropic, False if the problem is anisotropic
homogeneous = False # True if the growth rates are homogeneous, False if the problem is heterogeneous
useFibres = True # True if fibres are used for anisotropic problems
heterogeneousFibres = False # True if the fibre angles vary in space, False if not.
fixBottomRing = True # True if the bottom ring of nodes is fixed, False if not
fixTopRing = False # True if the top ring of nodes is fixed, False if not

# Tube geometry
length = 6.0 # The length of the tube
innerRadius = 1.25 # The inner radius of the tube
outerRadius = 2.0 # The outer radius of the tube

numberOfLengthElements=4 # Number of elements along the length of the tube
numberOfCircumfrentialElementsPerQuarter=1 # Number of elements in the circumfrential direction in one quarter of the tube
numberOfWallElements=1 # Number of elements through the wall of the tube

# Growth rates
baseGrowthRate = 0.005 # The base (isotropic) growth rate

fibreGrowthRatio = 1.0 # The multiplicative factor in the fibre direction applied to the base growth rate
sheetGrowthRatio = 1.0 # The multiplicative factor in the sheet direction applied to the base growth rate
normalGrowthRatio = 0.0 # The multiplicative factor in the normal direction applied to the base growth rate

heterogeneousFibreFactor = 1.0 # For fibre heteregeous problems, the fibre growth rates are multiplied by one plus the factor multiplied by the length fraction
heterogeneousSheetFactor = 3.0 # For heteregeous problems, the sheet growth rates are multiplied one plus the factor multiplied by cos(theta)
heterogeneousNormalFactor = 0.0 # For heteregeous problems, the normal growth rates are multiplied by one plus the factor multiplied by the wall fraction

# Hydrostatic pressure
pInit = -8.0 # The initial hydrostatic pressure

# Fibre angle
#fibreAngle = math.pi/2.0 # The fibre angle wrt the for anisotropic fibres
fibreAngle = 0.0 # The fibre angle wrt the for anisotropic fibres


# Fitting smoothing parameters
tau = 0.1 # Tension smoothing parameter for fitting
kappa = 0.05 # Curvature smoothing parameter for fitting

# Times
startTime = 0.0 # The start time for the growth simulation
stopTime1 = 30.0 # The stop time for the growth simulation
stopTime2 = 90.0 # The stop time for the growth simulation
timeIncrement = 1.0 # The time increment for the growth simulation

# Number of Gauss points used
numberOfGaussXi = 3

# Should not need to change anything below here

def SetGrowthPatternParameters(growthPhase):
    # Initialise the parameters field
    if not homogeneous:
        for wallElementIdx in range(1,numberOfWallElements+1):
            for lengthElementIdx in range(1,numberOfLengthElements+1):
                for circumfrentialElementIdx in range(1,numberOfCircumfrentialElements+1):
                    elementNumber = circumfrentialElementIdx + (lengthElementIdx-1)*numberOfCircumfrentialElements + \
                                    (wallElementIdx-1)*numberOfCircumfrentialElements*numberOfLengthElements
                    for xiIdx3 in range(1,numberOfGaussXi+1):
                        for xiIdx2 in range(1,numberOfGaussXi+1):
                            for xiIdx1 in range(1,numberOfGaussXi+1):
                                gaussPointNumber = xiIdx1 + (xiIdx2-1)*numberOfGaussXi + (xiIdx3-1)*numberOfGaussXi*numberOfGaussXi
                                radius = (float(wallElementIdx-1)+float(xiIdx3)/float(numberOfGaussXi+1))/ \
                                         float(numberOfWallElements)
                                theta = (float(circumfrentialElementIdx-1)+float(xiIdx1)/float(numberOfGaussXi+1))/ \
                                        float(numberOfCircumfrentialElements)*2.0*math.pi    
                                length = (float(lengthElementIdx-1)+float(xiIdx2)/float(numberOfGaussXi+1))/ \
                                         float(numberOfLengthElements)
                                if growthPhase == 1:
                                    fibreRate = fibreGrowthRate*heterogeneousFibreFactor
                                    sheetRate = sheetGrowthRate*heterogeneousSheetFactor*(1.0 + math.cos(theta))
                                    normalRate = normalGrowthRate*heterogeneousNormalFactor*(1.0 + radius)
                                elif growthPhase == 2:
                                    fibreRate = fibreGrowthRate*heterogeneousFibreFactor
                                    sheetRate = sheetGrowthRate*heterogeneousSheetFactor*(1.0 + math.sin(theta)*length*length)
                                    normalRate = normalGrowthRate*heterogeneousNormalFactor*(1.0 + radius)
                                else:
                                    print("The specified growth phase of "+growthPhase+" is invalid.")
                                growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                                                                           iron.FieldParameterSetTypes.VALUES,
                                                                                           gaussPointNumber,elementNumber,1,
                                                                                           fibreRate)
                                growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                                                                           iron.FieldParameterSetTypes.VALUES,
                                                                                           gaussPointNumber,elementNumber,2,
                                                                                           sheetRate)
                                growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                                                                           iron.FieldParameterSetTypes.VALUES,
                                                                                           gaussPointNumber,elementNumber,3,
                                                                                           normalRate)
                                growthCellMLParametersField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                                                                    iron.FieldParameterSetTypes.VALUES)
                                growthCellMLParametersField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                                                                     iron.FieldParameterSetTypes.VALUES)


    else:
        growthCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                                iron.FieldParameterSetTypes.VALUES,1,fibreGrowthRate)
        growthCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                                iron.FieldParameterSetTypes.VALUES,2,sheetGrowthRate)
        growthCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                                iron.FieldParameterSetTypes.VALUES,3,normalGrowthRate)


def SolveGrowthAndFittingProblems(time):
    timeLoop.TimesSet(time,time+timeIncrement,timeIncrement)
    # Solve the problem
    try:
        problem.Solve()
    except:
        print("Error during elasticity solve.\n")
    #Update lambda
    for wallElementIdx in range(1,numberOfWallElements+1):
        for lengthElementIdx in range(1,numberOfLengthElements+1):
            for circumfrentialElementIdx in range(1,numberOfCircumfrentialElements+1):
                elementNumber = circumfrentialElementIdx + (lengthElementIdx-1)*numberOfCircumfrentialElements + \
                                (wallElementIdx-1)*numberOfCircumfrentialElements*numberOfLengthElements
                for xiIdx3 in range(1,numberOfGaussXi+1):
                    for xiIdx2 in range(1,numberOfGaussXi+1):
                        for xiIdx1 in range(1,numberOfGaussXi+1):
                            gaussPointNumber = xiIdx1 + (xiIdx2-1)*numberOfGaussXi + (xiIdx3-1)*numberOfGaussXi*numberOfGaussXi
			    lambda1 = growthCellMLStateField.ParameterSetGetGaussPointDP(iron.FieldVariableTypes.U,
				iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,1)-1.0
			    fittingIndependentField.ParameterSetAddGaussPointDP(iron.FieldVariableTypes.U,
			    	iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,1,lambda1)
			    lambda2 = growthCellMLStateField.ParameterSetGetGaussPointDP(iron.FieldVariableTypes.U,
			    	iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,2)-1.0
			    fittingIndependentField.ParameterSetAddGaussPointDP(iron.FieldVariableTypes.U,
			    	iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,2,lambda2)
    # Update the independent field
    fittingIndependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    fittingIndependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    # Reset the fitting problem boundary conditions
    x = dependentField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,1,1)
    y = dependentField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,1,1)
    newR = math.sqrt(x*x+y*y)
    fittingDependentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                   1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,1,1,newR/innerRadius)
    fittingDependentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                   1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,1,2,1.0)
    # Update the independent field
    fittingDependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    fittingDependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    # Solve the fitting problem
    try:
        fittingProblem.Solve()
    except:
        print("Error during fitting solve.\n")
    #Copy fitted dependent field to lambda field
    iron.Field.ParametersToFieldParametersComponentCopy(
        fittingDependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,
        lambdaField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1)
    iron.Field.ParametersToFieldParametersComponentCopy(
        fittingDependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,
        lambdaField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2)
    # Export results
    timeString = format(time)
    filename = "HeartTubeGrowth_"+timeString
    fields.NodesExport(filename,"FORTRAN")
    fields.ElementsExport(filename,"FORTRAN")
    # Set geometric field to current deformed geometry
    iron.Field.ParametersToFieldParametersComponentCopy(
        dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,
        geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1)
    iron.Field.ParametersToFieldParametersComponentCopy(
        dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,
        geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2)
    iron.Field.ParametersToFieldParametersComponentCopy(
        dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,
        geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3)
    # Reset growth state to 1.0
    growthCellMLStateField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                       iron.FieldParameterSetTypes.VALUES,1,1.0)
    growthCellMLStateField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                       iron.FieldParameterSetTypes.VALUES,2,1.0)
    growthCellMLStateField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                       iron.FieldParameterSetTypes.VALUES,3,1.0)

if isotropic == True:
    fibreGrowthRate = baseGrowthRate
    sheetGrowthRate = baseGrowthRate
    normalGrowthRate = baseGrowthRate
else:
    fibreGrowthRate = fibreGrowthRatio*baseGrowthRate
    sheetGrowthRate = sheetGrowthRatio*baseGrowthRate
    normalGrowthRate = normalGrowthRatio*baseGrowthRate

numberOfCircumfrentialElements = 4*numberOfCircumfrentialElementsPerQuarter
numberOfLengthNodes = numberOfLengthElements+1
numberOfCircumfrentialNodes = numberOfCircumfrentialElements
numberOfWallNodes = numberOfWallElements+1

coordinateSystemUserNumber = 1
regionUserNumber = 1
tricubicHermiteBasisUserNumber = 1
trilinearLagrangeBasisUserNumber = 2
meshUserNumber = 1
decompositionUserNumber = 1
geometricFieldUserNumber = 1
fibreFieldUserNumber = 2
dependentFieldUserNumber = 3
equationsSetUserNumber = 1
equationsSetFieldUserNumber = 5
growthCellMLUserNumber = 1
growthCellMLModelsFieldUserNumber = 6
growthCellMLStateFieldUserNumber = 7
growthCellMLParametersFieldUserNumber = 8
constitutiveCellMLUserNumber = 2
constitutiveCellMLModelsFieldUserNumber = 9
constitutiveCellMLParametersFieldUserNumber = 10
constitutiveCellMLIntermediateFieldUserNumber = 11
problemUserNumber = 1
lambdaFieldUserNumber = 12
fittingEquationsSetUserNumber = 13
fittingEquationsSetFieldUserNumber = 14
fittingDependentFieldUserNumber = 15
fittingIndependentFieldUserNumber = 16
fittingMaterialsFieldUserNumber = 17
fittingProblemUserNumber = 18

#iron.DiagnosticsSetOn(iron.DiagnosticTypes.FROM,[1,2,3,4,5],"diagnostics",["FiniteElasticity_FiniteElementResidualEvaluate"])

# Get the number of computational nodes and this computational node number
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

# Create a 3D rectangular cartesian coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
# Set the number of dimensions to 3
coordinateSystem.DimensionSet(3)
# Finish the creation of the coordinate system
coordinateSystem.CreateFinish()

# Create a region and assign the coordinate system to the region
region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.LabelSet("HeartTubeRegion")
# Set the regions coordinate system to the 3D RC coordinate system that we have created
region.coordinateSystem = coordinateSystem
# Finish the creation of the region
region.CreateFinish()

# Define basis
# Start the creation of a tricubic Hermite basis function
tricubicHermiteBasis = iron.Basis()
tricubicHermiteBasis.CreateStart(tricubicHermiteBasisUserNumber)
tricubicHermiteBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
tricubicHermiteBasis.numberOfXi = 3
tricubicHermiteBasis.interpolationXi = [iron.BasisInterpolationSpecifications.CUBIC_HERMITE]*3
tricubicHermiteBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
tricubicHermiteBasis.CreateFinish()
# Start the creation of a trilinear Hermite basis function
trilinearLagrangeBasis = iron.Basis()
trilinearLagrangeBasis.CreateStart(trilinearLagrangeBasisUserNumber)
trilinearLagrangeBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
trilinearLagrangeBasis.numberOfXi = 3
trilinearLagrangeBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
trilinearLagrangeBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
trilinearLagrangeBasis.CreateFinish()

# Start the creation of a manually generated mesh in the region
numberOfNodes = numberOfCircumfrentialElements*(numberOfLengthElements+1)*(numberOfWallElements+1)
numberOfElements = numberOfCircumfrentialElements*numberOfLengthElements*numberOfWallElements


# Define nodes for the mesh
nodes = iron.Nodes()
nodes.CreateStart(region,numberOfNodes)
nodes.CreateFinish()
mesh = iron.Mesh()

# Create the mesh. The mesh will have two components - 1. tricubic Hermite elements; 2. trilinear Lagrange elements
mesh.CreateStart(meshUserNumber,region,3)
mesh.NumberOfComponentsSet(2)
mesh.NumberOfElementsSet(numberOfElements)

tricubicHermiteElements = iron.MeshElements()
tricubicHermiteElements.CreateStart(mesh,1,tricubicHermiteBasis)
trilinearLagrangeElements = iron.MeshElements()
trilinearLagrangeElements.CreateStart(mesh,2,trilinearLagrangeBasis)

elementNumber = 0
for wallElementIdx in range(1,numberOfWallElements+1):
    for lengthElementIdx in range(1,numberOfLengthElements+1):
        for circumfrentialElementIdx in range(1,numberOfCircumfrentialElements+1):
            elementNumber = elementNumber + 1
            localNode1 = circumfrentialElementIdx + (lengthElementIdx-1)*numberOfCircumfrentialNodes + \
                (wallElementIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
            if circumfrentialElementIdx == numberOfCircumfrentialElements:
                localNode2 = 1 + (lengthElementIdx-1)*numberOfCircumfrentialNodes + \
                    (wallElementIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
            else:
                localNode2 = localNode1 + 1
            localNode3 = localNode1 + numberOfCircumfrentialNodes
            localNode4 = localNode2 + numberOfCircumfrentialNodes
            localNode5 = localNode1 + numberOfCircumfrentialNodes*numberOfLengthNodes
            localNode6 = localNode2 + numberOfCircumfrentialNodes*numberOfLengthNodes
            localNode7 = localNode3 + numberOfCircumfrentialNodes*numberOfLengthNodes
            localNode8 = localNode4 + numberOfCircumfrentialNodes*numberOfLengthNodes
            localNodes = [localNode1,localNode2,localNode3,localNode4,localNode5,localNode6,localNode7,localNode8]
            tricubicHermiteElements.NodesSet(elementNumber,localNodes)
            trilinearLagrangeElements.NodesSet(elementNumber,localNodes)

tricubicHermiteElements.CreateFinish()
trilinearLagrangeElements.CreateFinish()

# Finish the mesh creation
mesh.CreateFinish() 

# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
# Set the decomposition to be a general decomposition with the specified number of domains
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
# Finish the decomposition
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
# Set the decomposition to use
geometricField.MeshDecompositionSet(decomposition)
geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
# Set the field label
geometricField.VariableLabelSet(iron.FieldVariableTypes.U,"Geometry")
# Set the domain to be used by the field components to be tricubic Hermite
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
# Set the scaling type
geometricField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
# Finish creating the field
geometricField.CreateFinish()

# Create the geometric field
for wallNodeIdx in range(1,numberOfWallNodes+1):
    for lengthNodeIdx in range(1,numberOfLengthNodes+1):
        for circumfrentialNodeIdx in range(1,numberOfCircumfrentialNodes+1):
            nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*numberOfCircumfrentialNodes + \
                (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
            radius = innerRadius + (outerRadius - innerRadius)*float(wallNodeIdx-1)/float(numberOfWallNodes)
            theta = float(circumfrentialNodeIdx-1)/float(numberOfCircumfrentialNodes)*2.0*math.pi
            x = radius*math.cos(theta)
            y = radius*math.sin(theta)
            xtangent = -math.sin(theta)
            ytangent = math.cos(theta)
            xnormal = math.cos(theta)
            ynormal = math.sin(theta)
            z = float(lengthNodeIdx-1)/float(numberOfLengthNodes)*length
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,x)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,y)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,z)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,xtangent)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,ytangent)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,3,0.0)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,0.0)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,0.0)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,3,1.0)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,1,xnormal)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,2,ynormal)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,3,0.0)

# Update the geometric field
geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

if useFibres:
    # Create a fibre field and attach it to the geometric field
    fibreField = iron.Field()
    fibreField.CreateStart(fibreFieldUserNumber,region)
    fibreField.TypeSet(iron.FieldTypes.FIBRE)
    # Set the decomposition 
    fibreField.MeshDecompositionSet(decomposition)
    # Set the geometric field
    fibreField.GeometricFieldSet(geometricField)
    # Set the field variable label
    fibreField.VariableLabelSet(iron.FieldVariableTypes.U,"Fibre")
    # Set the fibre field to use trilinear-Lagrange elements
    fibreField.NumberOfComponentsSet(iron.FieldVariableTypes.U,3)
    fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,2)
    fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,2)
    fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,2)
    # Finish creating the field
    fibreField.CreateFinish()
    #Initialise the fibre field
    for wallNodeIdx in range(1,numberOfWallNodes+1):
        for lengthNodeIdx in range(1,numberOfLengthNodes+1):
            for circumfrentialNodeIdx in range(1,numberOfCircumfrentialNodes+1):
                nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*numberOfCircumfrentialNodes + \
                    (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
                # Set the fibre angle
                if heterogeneousFibres == True:
            	    theta = float(circumfrentialNodeIdx-1)/float(numberOfCircumfrentialNodes)*2.0*math.pi
		    angle = fibreAngle*math.sin(theta)
                else:
		    angle = fibreAngle
                fibreField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                        1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,angle)
                fibreField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                        1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,0.0)
                fibreField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                        1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,0.0)
    # Update the fibre field
    fibreField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    fibreField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

# Export results
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("HeartTubeGrowth","FORTRAN")
fields.ElementsExport("HeartTubeGrowth","FORTRAN")
fields.Finalise()

# Create the dependent field
dependentField = iron.Field()
dependentField.CreateStart(dependentFieldUserNumber,region)
dependentField.TypeSet(iron.FieldTypes.GEOMETRIC_GENERAL)
# Set the decomposition
dependentField.MeshDecompositionSet(decomposition)
# Set the geometric field
dependentField.GeometricFieldSet(geometricField) 
dependentField.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT)
# Set the field variables for displacement, traction, strain, stress and growth
dependentField.NumberOfVariablesSet(5)
dependentField.VariableTypesSet([iron.FieldVariableTypes.U,iron.FieldVariableTypes.DELUDELN,
                                 iron.FieldVariableTypes.U1,iron.FieldVariableTypes.U2,iron.FieldVariableTypes.U3])
dependentField.VariableLabelSet(iron.FieldVariableTypes.U,"Dependent")
dependentField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"del U/del n")
dependentField.VariableLabelSet(iron.FieldVariableTypes.U1,"Strain")
dependentField.VariableLabelSet(iron.FieldVariableTypes.U2,"Stress")
dependentField.VariableLabelSet(iron.FieldVariableTypes.U3,"Growth")
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U,4)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN,4)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U1,6)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U2,6)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U3,3)
# Set the hydrostatic pressure to use tri-linear Lagrange elements
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,4,2)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,4,2)
# Set the strain, stress and growth variables to be Gauss point based.
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,1,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,2,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,3,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,4,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,5,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,6,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,1,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,2,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,3,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,4,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,5,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,6,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U3,1,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U3,2,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U3,3,
                                         iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
# Set the field scaling
dependentField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
# Finish creating the field
dependentField.CreateFinish()

# Initialise dependent field from undeformed geometry
iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,
    dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1)
iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,
    dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2)
iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,
    dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3)
# Initialise the hydrostatic pressure
iron.Field.ComponentValuesInitialiseDP(dependentField,iron.FieldVariableTypes.U,
                                       iron.FieldParameterSetTypes.VALUES,4,pInit)

# Create the equations_set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
# Specify a finite elasticity equations set with the growth and constitutive law in CellML
equationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
    iron.EquationsSetTypes.FINITE_ELASTICITY,
    iron.EquationsSetSubtypes.CONSTIT_AND_GROWTH_LAW_IN_CELLML]
if useFibres:
    equationsSet.CreateStart(equationsSetUserNumber,region,fibreField,
                             equationsSetSpecification,equationsSetFieldUserNumber,
                             equationsSetField)
else:
    equationsSet.CreateStart(equationsSetUserNumber,region,geometricField,
                             equationsSetSpecification,equationsSetFieldUserNumber,
                             equationsSetField)
equationsSet.CreateFinish()

# Set up the equation set dependent field
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
equationsSet.DependentCreateFinish()

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
# Use sparse equations
equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
# Do not output any equations information
equations.outputType = iron.EquationsOutputTypes.NONE
# Finish creating the equations
equationsSet.EquationsCreateFinish()

# Set up the growth CellML model
growthCellML = iron.CellML()
growthCellML.CreateStart(growthCellMLUserNumber,region)
if growthModel == 1:
   # Create the CellML environment for the simple growth law
   growthCellMLIdx = growthCellML.ModelImport("simplegrowth.cellml")
   # Flag the CellML variables that OpenCMISS will supply
   growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/fibrerate")
   growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/sheetrate")
   growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/normalrate")
elif growthModel == 2:
   # Create the CellML environment for the stress based growth law
   growthCellML = iron.CellML()
   growthCellML.CreateStart(growthCellMLUserNumber,region)
   growthCellMLIdx = growthCellML.ModelImport("stressgrowth.cellml")
   # Flag the CellML variables that OpenCMISS will supply
   growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/fibrerate")
   growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/sheetrate")
   growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/normalrate")
   growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/S11")
   growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/S22")
   growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/S33")
elif growthModel == 3:
   # Create the CellML environment for the strain based growth law
   growthCellML = iron.CellML()
   growthCellML.CreateStart(growthCellMLUserNumber,region)
   growthCellMLIdx = growthCellML.ModelImport("straingrowth.cellml")
   # Flag the CellML variables that OpenCMISS will supply
   growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/fibrerate")
   growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/sheetrate")
   growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/normalrate")
   growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/C11")
   growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/C22")
   growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/C33")
# Finish the growth CellML
growthCellML.CreateFinish()

# Create CellML <--> OpenCMISS field maps
growthCellML.FieldMapsCreateStart()
if growthModel == 2:
   growthCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U2,1,iron.FieldParameterSetTypes.VALUES,
	growthCellMLIdx,"Main/S11",iron.FieldParameterSetTypes.VALUES)
   growthCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U2,2,iron.FieldParameterSetTypes.VALUES,
	growthCellMLIdx,"Main/S22",iron.FieldParameterSetTypes.VALUES)
   growthCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U3,3,iron.FieldParameterSetTypes.VALUES,
	growthCellMLIdx,"Main/S33",iron.FieldParameterSetTypes.VALUES)
elif growthModel == 3:
   growthCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,1,iron.FieldParameterSetTypes.VALUES,
	growthCellMLIdx,"Main/C11",iron.FieldParameterSetTypes.VALUES)
   growthCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,2,iron.FieldParameterSetTypes.VALUES,
	growthCellMLIdx,"Main/C22",iron.FieldParameterSetTypes.VALUES)
   growthCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,3,iron.FieldParameterSetTypes.VALUES,
	growthCellMLIdx,"Main/C33",iron.FieldParameterSetTypes.VALUES)
growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda1",iron.FieldParameterSetTypes.VALUES,
	dependentField,iron.FieldVariableTypes.U3,1,iron.FieldParameterSetTypes.VALUES)
growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda2",iron.FieldParameterSetTypes.VALUES,
	dependentField,iron.FieldVariableTypes.U3,2,iron.FieldParameterSetTypes.VALUES)
growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda3",iron.FieldParameterSetTypes.VALUES,
        dependentField,iron.FieldVariableTypes.U3,3,iron.FieldParameterSetTypes.VALUES)
growthCellML.FieldMapsCreateFinish()

# Create the CELL models field
growthCellMLModelsField = iron.Field()
growthCellML.ModelsFieldCreateStart(growthCellMLModelsFieldUserNumber,growthCellMLModelsField)
growthCellMLModelsField.VariableLabelSet(iron.FieldVariableTypes.U,"GrowthModelMap")
growthCellML.ModelsFieldCreateFinish()

# Create the CELL parameters field
growthCellMLParametersField = iron.Field()
growthCellML.ParametersFieldCreateStart(growthCellMLParametersFieldUserNumber,growthCellMLParametersField)
growthCellMLParametersField.VariableLabelSet(iron.FieldVariableTypes.U,"GrowthParameters")
growthCellML.ParametersFieldCreateFinish()


# Create the CELL state field
growthCellMLStateField = iron.Field()
growthCellML.StateFieldCreateStart(growthCellMLStateFieldUserNumber,growthCellMLStateField)
growthCellMLStateField.VariableLabelSet(iron.FieldVariableTypes.U,"GrowthState")
growthCellML.StateFieldCreateFinish()

# Create the CellML environment for the consitutative law
constitutiveCellML = iron.CellML()
constitutiveCellML.CreateStart(constitutiveCellMLUserNumber,region)
constitutiveCellMLIdx = constitutiveCellML.ModelImport("mooneyrivlin.cellml")
# Flag the CellML variables that OpenCMISS will supply
constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/E11")
constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/E12")
constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/E13")
constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/E22")
constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/E23")
constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/E33")
#constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/c1")
#constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/c2")
# Flag the CellML variables that OpenCMISS will obtain
constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev11")
constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev12")
constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev13")
constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev22")
constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev23")
constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev33")
constitutiveCellML.CreateFinish()

# Create CellML <--> OpenCMISS field maps
constitutiveCellML.FieldMapsCreateStart()
constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,1,iron.FieldParameterSetTypes.VALUES,
    constitutiveCellMLIdx,"equations/E11",iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,2,iron.FieldParameterSetTypes.VALUES,
    constitutiveCellMLIdx,"equations/E12",iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,3,iron.FieldParameterSetTypes.VALUES,
    constitutiveCellMLIdx,"equations/E13",iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,4,iron.FieldParameterSetTypes.VALUES,
    constitutiveCellMLIdx,"equations/E22",iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,5,iron.FieldParameterSetTypes.VALUES,
    constitutiveCellMLIdx,"equations/E23",iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,6,iron.FieldParameterSetTypes.VALUES,
    constitutiveCellMLIdx,"equations/E33",iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev11",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U2,1,iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev12",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U2,2,iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev13",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U2,3,iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev22",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U2,4,iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev23",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U2,5,iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev33",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U2,6,iron.FieldParameterSetTypes.VALUES)
constitutiveCellML.FieldMapsCreateFinish()

# Create the CELL models field
constitutiveCellMLModelsField = iron.Field()
constitutiveCellML.ModelsFieldCreateStart(constitutiveCellMLModelsFieldUserNumber,
                                           constitutiveCellMLModelsField)
constitutiveCellMLModelsField.VariableLabelSet(iron.FieldVariableTypes.U,"ConstitutiveModelMap")
constitutiveCellML.ModelsFieldCreateFinish()

# Create the CELL parameters field
constitutiveCellMLParametersField = iron.Field()
constitutiveCellML.ParametersFieldCreateStart(constitutiveCellMLParametersFieldUserNumber,
                                               constitutiveCellMLParametersField)
constitutiveCellMLParametersField.VariableLabelSet(iron.FieldVariableTypes.U,"ConstitutiveParameters")
constitutiveCellML.ParametersFieldCreateFinish()

# Create the CELL intermediate field
constitutiveCellMLIntermediateField = iron.Field()
constitutiveCellML.IntermediateFieldCreateStart(constitutiveCellMLIntermediateFieldUserNumber,
                                                 constitutiveCellMLIntermediateField)
constitutiveCellMLIntermediateField.VariableLabelSet(iron.FieldVariableTypes.U,"ConstitutiveIntermediate")
constitutiveCellML.IntermediateFieldCreateFinish()

# Define the problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.ELASTICITY,
        iron.ProblemTypes.FINITE_ELASTICITY,
        iron.ProblemSubtypes.FINITE_ELASTICITY_WITH_GROWTH_CELLML]
problem.CreateStart(problemUserNumber,problemSpecification)
problem.CreateFinish()

# Create control loops
timeLoop = iron.ControlLoop()
problem.ControlLoopCreateStart()
problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],timeLoop)
problem.ControlLoopCreateFinish()

# Create problem solvers
odeIntegrationSolver = iron.Solver()
nonlinearSolver = iron.Solver()
linearSolver = iron.Solver()
cellMLEvaluationSolver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,odeIntegrationSolver)
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,nonlinearSolver)
#nonlinearSolver.outputType = iron.SolverOutputTypes.MONITOR
nonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD)
nonlinearSolver.NewtonCellMLSolverGet(cellMLEvaluationSolver)
nonlinearSolver.NewtonLinearSolverGet(linearSolver)
linearSolver.linearType = iron.LinearSolverTypes.DIRECT
problem.SolversCreateFinish()

# Create nonlinear equations and add equations set to solver equations
nonlinearEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
nonlinearSolver.SolverEquationsGet(nonlinearEquations)
nonlinearEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
nonlinearEquationsSetIndex = nonlinearEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Create CellML equations and add growth and constitutive equations to the solvers
growthEquations = iron.CellMLEquations()
constitutiveEquations = iron.CellMLEquations()
problem.CellMLEquationsCreateStart()
odeIntegrationSolver.CellMLEquationsGet(growthEquations)
growthEquationsIndex = growthEquations.CellMLAdd(growthCellML)
cellMLEvaluationSolver.CellMLEquationsGet(constitutiveEquations)
constitutiveEquationsIndex = constitutiveEquations.CellMLAdd(constitutiveCellML)
problem.CellMLEquationsCreateFinish()

# Prescribe boundary conditions (absolute nodal parameters)
boundaryConditions = iron.BoundaryConditions()
nonlinearEquations.BoundaryConditionsCreateStart(boundaryConditions)

for lengthNodeIdx in range(1,3):
    if (lengthNodeIdx == 1 and fixBottomRing) or (lengthNodeIdx == 2 and fixTopRing):
        for wallNodeIdx in range(1,numberOfWallNodes+1):
            for circumfrentialNodeIdx in range(1,numberOfCircumfrentialNodes+1):
                nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*(numberOfLengthNodes-1)*numberOfCircumfrentialNodes + \
                             (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
                # Fix z direction
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                           1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                # Fix S1 (circumfrential) direction derivatives
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                           1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                           1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                           1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,3,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                # Fix S2 (length) direction derivatives
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                           1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                           1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                           1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,3,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                # Fix S3 (wall) direction derivatives
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                           1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,1,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                           1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,2,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                           1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,3,
                                           iron.BoundaryConditionsTypes.FIXED,0.0)
        #Set symmetry conditions on the ring to prevent rotation
        nodeNumber = 1 + (lengthNodeIdx-1)*(numberOfLengthNodes-1)*numberOfCircumfrentialNodes + \
                     (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                   1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,
                                   iron.BoundaryConditionsTypes.FIXED,0.0)
        nodeNumber = nodeNumber + numberOfCircumfrentialElementsPerQuarter
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                   1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,
                                   iron.BoundaryConditionsTypes.FIXED,0.0)
        nodeNumber = nodeNumber + numberOfCircumfrentialElementsPerQuarter
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                   1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,
                                   iron.BoundaryConditionsTypes.FIXED,0.0)
        nodeNumber = nodeNumber + numberOfCircumfrentialElementsPerQuarter
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                   1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,
                                   iron.BoundaryConditionsTypes.FIXED,0.0)

nonlinearEquations.BoundaryConditionsCreateFinish()

# Create the lambda field
lambdaField = iron.Field()
lambdaField.CreateStart(lambdaFieldUserNumber,region)
lambdaField.TypeSet(iron.FieldTypes.GENERAL)
# Set the decomposition
lambdaField.MeshDecompositionSet(decomposition)
# Set the geometric field
lambdaField.GeometricFieldSet(geometricField)
lambdaField.ScalingTypeSet(iron.FieldScalingTypes.NONE)
# Set the field variables
lambdaField.NumberOfVariablesSet(1)
lambdaField.VariableTypesSet([iron.FieldVariableTypes.U])
# Set the variable label
lambdaField.VariableLabelSet(iron.FieldVariableTypes.U,"NodeLambda")
# Set the components to be trilinear-Lagrange
lambdaField.NumberOfComponentsSet(iron.FieldVariableTypes.U,2)
lambdaField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,2)
lambdaField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,2)
# Set the interpolation types
lambdaField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
lambdaField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
lambdaField.CreateFinish()
# Initialise the lambda field
lambdaField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1.0)
lambdaField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,1.0)

# Create Gauss point fitting equations set
fittingEquationsSetSpecification = [iron.EquationsSetClasses.FITTING,
                             iron.EquationsSetTypes.GAUSS_FITTING_EQUATION,
                             iron.EquationsSetSubtypes.GAUSS_POINT_FITTING,
                             iron.EquationsSetFittingSmoothingTypes.SOBOLEV_VALUE]
fittingEquationsSetField = iron.Field()
fittingEquationsSet = iron.EquationsSet()
fittingEquationsSet.CreateStart(fittingEquationsSetUserNumber,region,geometricField,
        fittingEquationsSetSpecification,fittingEquationsSetFieldUserNumber,fittingEquationsSetField)
fittingEquationsSet.CreateFinish()

# Create the fitting dependent field
fittingDependentField = iron.Field()
fittingEquationsSet.DependentCreateStart(fittingDependentFieldUserNumber,fittingDependentField)
fittingDependentField.VariableLabelSet(iron.FieldVariableTypes.U,"FittingU")
fittingDependentField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"FittingDelUdelN")
# Set the number of components to 2
fittingDependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U,2)
fittingDependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN,2)
# Set the field variables to be trilinear Lagrange
fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,2)
fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,2)
fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,1,2)
fittingDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,2,2)
# Finish creating the fitting dependent field
fittingEquationsSet.DependentCreateFinish()

# Create the fitting independent field
fittingIndependentField = iron.Field()
fittingEquationsSet.IndependentCreateStart(fittingIndependentFieldUserNumber,fittingIndependentField)
fittingIndependentField.VariableLabelSet(iron.FieldVariableTypes.U,"GaussLambda")
fittingIndependentField.VariableLabelSet(iron.FieldVariableTypes.V,"LambdaWeight")
# Set the number of components to 2
fittingIndependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U,2)
fittingIndependentField.NumberOfComponentsSet(iron.FieldVariableTypes.V,2)
# Finish creating the fitting independent field
fittingEquationsSet.IndependentCreateFinish()
# Initialise data point vector field to 0.0
fittingIndependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)
fittingIndependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,0.0)
# Initialise data point weight field to 1.0
fittingIndependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,1,1.0)
fittingIndependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,2,1.0)

# Create material field (Sobolev parameters)
fittingMaterialField = iron.Field()
fittingEquationsSet.MaterialsCreateStart(fittingMaterialsFieldUserNumber,fittingMaterialField)
fittingMaterialField.VariableLabelSet(iron.FieldVariableTypes.U,"SmoothingParameters")
fittingEquationsSet.MaterialsCreateFinish()
# Set kappa and tau - Sobolev smoothing parameters
fittingMaterialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,tau)
fittingMaterialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,kappa)

# Create the fitting equations
fittingEquations = iron.Equations()
fittingEquationsSet.EquationsCreateStart(fittingEquations)
# Set the fitting equations sparsity type
fittingEquations.sparsityType = iron.EquationsSparsityTypes.SPARSE
# Set the fitting equations output type to none
fittingEquations.outputType = iron.EquationsOutputTypes.NONE
# Finish creating the fitting equations
fittingEquationsSet.EquationsCreateFinish()

# Create fitting problem
fittingProblemSpecification = [iron.ProblemClasses.FITTING,
                        iron.ProblemTypes.DATA_FITTING,
                        iron.ProblemSubtypes.STATIC_FITTING]
fittingProblem = iron.Problem()
fittingProblem.CreateStart(fittingProblemUserNumber,fittingProblemSpecification)
fittingProblem.CreateFinish()

# Create control loops
fittingProblem.ControlLoopCreateStart()
fittingProblem.ControlLoopCreateFinish()

# Create problem solver
fittingSolver = iron.Solver()
fittingProblem.SolversCreateStart()
fittingProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,fittingSolver)
fittingSolver.outputType = iron.SolverOutputTypes.PROGRESS
fittingProblem.SolversCreateFinish()

# Create fitting solver equations and add fitting equations set to solver equations
fittingSolverEquations = iron.SolverEquations()
fittingProblem.SolverEquationsCreateStart()
# Get the solver equations
fittingSolver.SolverEquationsGet(fittingSolverEquations)
fittingSolverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
fittingEquationsSetIndex = fittingSolverEquations.EquationsSetAdd(fittingEquationsSet)
fittingProblem.SolverEquationsCreateFinish()

# Prescribe boundary conditions for the fitting problem
fittingBoundaryConditions = iron.BoundaryConditions()
fittingSolverEquations.BoundaryConditionsCreateStart(fittingBoundaryConditions)
fittingBoundaryConditions.SetNode(fittingDependentField,iron.FieldVariableTypes.U,
	1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,1,1,
        iron.BoundaryConditionsTypes.FIXED,1.0)	
fittingSolverEquations.BoundaryConditionsCreateFinish()

SetGrowthPatternParameters(1)

# Loop over the time steps
time = startTime
timeString = format(time)
filename = "HeartTubeGrowth_"+timeString
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport(filename,"FORTRAN")
fields.ElementsExport(filename,"FORTRAN")

currentLength=0.0
while time<=stopTime1:
    SolveGrowthAndFittingProblems(time)
    time = time+timeIncrement
    
SetGrowthPatternParameters(2)

while time<=stopTime2:
    SolveGrowthAndFittingProblems(time)
    time = time+timeIncrement
    
fields.Finalise()

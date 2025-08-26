# Reference: https://opensees.berkeley.edu/wiki/index.php?title=Time_History_Analysis_of_a_2D_Elastic_Cantilever_Column
# --------------------------------------------------------------------------------------------------
# Example 1. cantilever 2D
# EQ ground motion with gravity
# all units are in kip, inch, second
# elasticBeamColumn ELEMENT
#		Silvia Mazzoni & Frank McKenna, 2006
# 		Marius Bittner, 2024

#SI-Units: seconds, meters, kilogramm, newtons
#
#    ^Y	
#    1   50m   |             
#    2         |
#    3         |  
#    4         | 
#    5         | 
#    6         | 
#    7         |
#    8         | 
#    9         | 
#    10        |
#    11        |   
#    12        |  
#    13        |
#    14        |
#    15        |
#    16        |
#    17  0m  --------->X


# SET UP ----------------------------------------------------------------------------
wipe;						       # clear opensees model
model basic -ndm 2 -ndf 3;	       #  2 = number of dimensions, 3 = degrees of freedom per node

# define GEOMETRY -------------------------------------------------------------
# nodal coordinates:
node 1  0.0 50.00
node 2  0.0 47.50
node 3  0.0 45.00
node 4  0.0 42.50
node 5  0.0 40.42
node 6  0.0 39.58
node 7  0.0 37.52
node 8  0.0 35.10
node 9  0.0 32.68
node 10 0.0 30.21
node 11 0.0 27.08
node 12 0.0 22.69
node 13 0.0 17.65
node 14 0.0 12.61
node 15 0.0  7.56
node 16 0.0  2.52
node 17 0.0  0.0

# Single point constraints -- Boundary Conditions
fix 17 1 1 1;  			           # node DX DY RZ (0 = unconstrained, 1 = constrained)

# nodal masses:
# Mass in the Y-direction is needed for gravity loads
mass 1 0 210.441   0.
mass 2 0 420.882   0.
mass 3 0 637.495   0.
mass 4 0 420.882   0.
mass 5 0 207.056   0.
mass 6 0 206.732   0.
mass 7 0 423.952   0.
mass 8 0 268.490   0.
mass 9 0 498.694   0.
mass 10 0 557.621  0.
mass 11 0 903.022  0.
mass 12 0 1359.248 0.
mass 13 0 1521.297 0.
mass 14 0 1683.347 0.
mass 15 0 1845.396 0.
mass 16 0 2007.446 0.               
                                       
# Define ELEMENTS -------------------------------------------------------------
# define geometric transformation: performs a linear geometric transformation of beam stiffness and resisting force from the basic system to the global-coordinate system
 # geomTransf Linear $transfTag     
geomTransf Linear 1;  
		      
#define section properties: section Elastic $secTag $E $A $Iz
section Elastic 1  {{{:E}}} 0.007 0.016
section Elastic 2  {{{:E}}} 0.007 0.016
section Elastic 3  {{{:E}}} 0.007 0.037
section Elastic 4  {{{:E}}} 0.007 0.016
section Elastic 5  {{{:E}}} 0.007 0.016
section Elastic 6  {{{:E}}} 0.011 0.023
section Elastic 7  {{{:E}}} 0.011 0.026
section Elastic 8  {{{:E}}} 0.011 0.011
section Elastic 9  {{{:E}}} 0.011 0.037
section Elastic 10 {{{:E}}} 0.014 0.054
section Elastic 11 {{{:E}}} 0.014 0.064
section Elastic 12 {{{:E}}} 0.014 0.079
section Elastic 13 {{{:E}}} 0.017 0.123
section Elastic 14 {{{:E}}} 0.017 0.150
section Elastic 15 {{{:E}}} 0.022 0.232
section Elastic 16 {{{:E}}} 0.022 0.275

# connectivity: element elasticBeamColumn $eleTag $iNode $jNode $secTag $transfTag
element elasticBeamColumn 1  1  2  1  1;
element elasticBeamColumn 2  2  3  2  1;
element elasticBeamColumn 3  3  4  3  1;
element elasticBeamColumn 4  4  5  4  1;
element elasticBeamColumn 5  5  6  5  1;
element elasticBeamColumn 6  6  7  6  1;
element elasticBeamColumn 7  7  8  7  1;
element elasticBeamColumn 8  8  9  8  1;
element elasticBeamColumn 9  9 10  9  1;
element elasticBeamColumn 10 10 11 10  1;
element elasticBeamColumn 11 11 12 11  1;
element elasticBeamColumn 12 12 13 12  1;
element elasticBeamColumn 13 13 14 13  1;
element elasticBeamColumn 14 14 15 14  1;
element elasticBeamColumn 15 15 16 15  1;
element elasticBeamColumn 16 16 17 16  1;

# define GRAVITY using UniformExcitation -------------------------------------
timeSeries Constant 1
# UniformExcitation pattern = apply a uniform excitation to a model acting in a certain direction
pattern UniformExcitation 1 2 -accel 1
#pattern UniformExcitation $patternTag $dir -accel $tsTag
# 1 = pattern tag
# 2 = direction (1=X, 2=Y)
# -accel -9.81 = gravity acceleration (negative for downward)

constraints Plain;          # how it handles boundary conditions    				
numberer Plain;	            # renumber dof's to minimize band-width (optimization), if you want to				    
system BandGeneral;		    # how to store and solve the system of equations in the analysi		    
algorithm Linear;           # use Linear algorithm for linear analysis              
integrator LoadControl 0.1;	# determine the next time step for an analysis		
analysis Static			    # define type of analysis: static		    
analyze 10;					# apply 10 load steps in analysis        
loadConst -time 0.0;        # set time to zero for static analysis

# Define RECORDERS -------------------------------------------------------------
# displacements of free nodes
# Node recorder type records the response of a number of nodes at every converged step
recorder Node -file displacement.out -time -node 1 -dof 1 disp		            
# recorder Node <-file $fileName><-time><-node $node1 $node2 ...>-dof ($dof1 $dof2 ...) $respType'
#disp = displacement

# DYNAMIC wind load analysis -------------------------------------------------------------
# define time series for wind load
#timeSeries = construct a TimeSeries object which represents the relationship between the time in the domain, and the load factor applied to the loads
#The input points can come from a file or from a list in the skript
#For a load path where the factors are specified in a file for a constant time interval between points:
#timeSeries Path $tag -dt $dt -filePath $filePath <-factor $cFactor>
#pattern Plain = construct a LoadPattern object $patternTag $tsTag <-fact $cFactor> {load}
#noad = command to nodal load  $nodeTag (ndf $LoadValues)
# Tags (2–17) are used for timeSeries and patterns to avoid conflicts
# Loads are applied in the X direction (1.0 0.0 0.0)
timeSeries Path 2  -dt {{{ :dt }}} -filePath wind-load-node1.dat  -factor 1.0;
pattern Plain 2 2 -timeSeries 2 {
    load 1 1.0 0.0 0.0;
}

timeSeries Path 3  -dt {{{ :dt }}} -filePath wind-load-node2.dat  -factor 1.0;
pattern Plain 3 3 -timeSeries 3 {
    load 2 1.0 0.0 0.0;
}

timeSeries Path 4  -dt {{{ :dt }}} -filePath wind-load-node3.dat  -factor 1.0;
pattern Plain 4 4 -timeSeries 4 {
    load 3 1.0 0.0 0.0;
}

timeSeries Path 5  -dt {{{ :dt }}} -filePath wind-load-node4.dat  -factor 1.0;
pattern Plain 5 5 -timeSeries 5 {
    load 4 1.0 0.0 0.0;
}

timeSeries Path 6  -dt {{{ :dt }}} -filePath wind-load-node5.dat  -factor 1.0;
pattern Plain 6 6 -timeSeries 6 {
    load 5 1.0 0.0 0.0;
}

timeSeries Path 7  -dt {{{ :dt }}} -filePath wind-load-node6.dat  -factor 1.0;
pattern Plain 7 7 -timeSeries 7 {
    load 6 1.0 0.0 0.0;
}

timeSeries Path 8  -dt {{{ :dt }}} -filePath wind-load-node7.dat  -factor 1.0;
pattern Plain 8 8 -timeSeries 8 {
    load 7 1.0 0.0 0.0;
}

timeSeries Path 9  -dt {{{ :dt }}} -filePath wind-load-node8.dat  -factor 1.0;
pattern Plain 9 9 -timeSeries 9 {
    load 8 1.0 0.0 0.0;
}

timeSeries Path 10 -dt {{{ :dt }}} -filePath wind-load-node9.dat  -factor 1.0;
pattern Plain 10 10 -timeSeries 10 {
    load 9 1.0 0.0 0.0;
}

timeSeries Path 11 -dt {{{ :dt }}} -filePath wind-load-node10.dat -factor 1.0;
pattern Plain 11 11 -timeSeries 11 {
    load 10 1.0 0.0 0.0;
}

timeSeries Path 12 -dt {{{ :dt }}} -filePath wind-load-node11.dat -factor 1.0;
pattern Plain 12 12 -timeSeries 12 {
    load 11 1.0 0.0 0.0;
}

timeSeries Path 13 -dt {{{ :dt }}} -filePath wind-load-node12.dat -factor 1.0;
pattern Plain 13 13 -timeSeries 13 {
    load 12 1.0 0.0 0.0;
}

timeSeries Path 14 -dt {{{ :dt }}} -filePath wind-load-node13.dat -factor 1.0;
pattern Plain 14 14 -timeSeries 14 {
    load 13 1.0 0.0 0.0;
}

timeSeries Path 15 -dt {{{ :dt }}} -filePath wind-load-node14.dat -factor 1.0;
pattern Plain 15 15 -timeSeries 15 {
    load 14 1.0 0.0 0.0;
}

timeSeries Path 16 -dt {{{ :dt }}} -filePath wind-load-node15.dat -factor 1.0;
pattern Plain 16 16 -timeSeries 16 {
    load 15 1.0 0.0 0.0;
}

timeSeries Path 17 -dt {{{ :dt }}} -filePath wind-load-node16.dat -factor 1.0;
pattern Plain 17 17 -timeSeries 17 {
    load 16 1.0 0.0 0.0;
}

# set damping based on first eigen mode
set freq [expr [eigen -fullGenLapack 1]**0.5]
set dampRatio 0.02
rayleigh 0. 0. 0. [expr 2*$dampRatio/$freq]

# create the analysis
wipeAnalysis;					     # clear previously-define analysis parameters
constraints Plain;     				 # how it handles boundary conditions
numberer Plain;					     # renumber dof's to minimize band-width (optimization), if you want to
system BandGeneral;					 # how to store and solve the system of equations in the analysis
algorithm Linear					 # use Linear algorithm for linear analysis
integrator Newmark 0.5 0.25 ;	     # determine the next time step for an analysis
analysis Transient;					 # define type of analysis: time-dependent
#analyze 3995 0.01;					 # apply 3995 0.01-sec time steps in analysis
# Read number of lines in wind-load-node2.dat (or node3)
set dt {{{ :dt }}}
set tEnd {{{ :T }}}
# Calculate the number of steps so that you get N = (tEnd/dt) + 1 output points
set nSteps [expr int($tEnd / $dt)]
record ; # Write initial state (t=0) to the output file
analyze $nSteps $dt

puts "Done!"
wipe
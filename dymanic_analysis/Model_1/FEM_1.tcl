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
#    |                  
#    |
#    3     50m _  
#    |         | 
#    |         | 
#    2     40m | 
#    |         |
#    |         | 
#    |         | 
#    |         | 
#    1    ----  -------->X
#

# SET UP ----------------------------------------------------------------------------
wipe;						       # clear opensees model
model basic -ndm 2 -ndf 3;	       #  2 = number of dimensions, 3 = degrees of freedom per node

# define GEOMETRY -------------------------------------------------------------
# nodal coordinates:
node 1 0. 0.;					   # node X Y
node 2 0. 40;                      # length in meters
node 3 0. 50;

# Single point constraints -- Boundary Conditions
fix 1 1 1 1; 			           # node DX DY RZ (0 = unconstrained, 1 = constrained)

# nodal masses:
mass 2 0 4.4e+03 0.;			   # node XD DY RZ (Mass=Weight in kg)
mass 3 0 8.8e+03 0.;               # Mass in the Y-direction is needed for gravity loads
                                       
# Define ELEMENTS -------------------------------------------------------------
# define geometric transformation: performs a linear geometric transformation of beam stiffness and resisting force from the basic system to the global-coordinate system
geomTransf Linear 1;  		       # geomTransf Linear $transfTag     

#define section properties: section Elastic $secTag $E $A $Iz
section Elastic 1 {{{:E}}} 0.1 {{{:Iz}}}  

# connectivity: element elasticBeamColumn $eleTag $iNode $jNode $secTag $transfTag
element elasticBeamColumn 1 1 2 1 1;  
element elasticBeamColumn 2 2 3 1 1;  


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
recorder Node -file displacement.out -time -node 3 -dof 1 disp		            
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

# For node 2
timeSeries Path 2 -dt {{{ :dt }}} -filePath wind-load-node2.dat -factor 1.0;
pattern Plain 2 2 -timeSeries 2 {
    load 2 1.0 0.0 0.0;
}

# For node 3
timeSeries Path 3 -dt {{{ :dt }}} -filePath wind-load-node3.dat -factor 1.0;
pattern Plain 3 3 -timeSeries 3 {
    load 3 1.0 0.0 0.0
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
set tEnd 40.0
# Calculate the number of steps so that you get N = (tEnd/dt) + 1 output points
set nSteps [expr int($tEnd / $dt)]
record ; # Write initial state (t=0) to the output file
analyze $nSteps $dt

puts "Done!"
wipe
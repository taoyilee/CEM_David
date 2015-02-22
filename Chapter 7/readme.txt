The files in this directory are MATLAB files, as discussed in Chapter 7,
D.B.Davidson, "Computational Electromagnetics for RF and Microwave Engineering",
cUP 2005, 2011. See the headers in the invididual files for details.

Main files are:

scalar_pot.m 
This script studies the behaviour of the scalar potential. 

V_pot_eps and V_pot_height
These scripts plot the behaviour of the scalar potential for varying relative permittivity and
height respecitively.

MoM_Som
This script computes the input impedance of a thin printed dipole over a frequency band, and permits users to specify both the number of modes and integration accuracy. It replaces MoM_Som_ver2.


The other files in the directory are functions required by the main scripts above. 

Known issues: Not all the codes work correctly with non-zero loss tangent, and this is not presently checked in the codes. For instance, the code generating Figs. 7.4 and 7.5 can handle lossy dieelectics, but to replicate Figs. 7.6-7.12, the loss should be set to zero. (The loss tangent is hard-wired in the code and must be edited).
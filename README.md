currently just solver for a bare square vacuum boundary, no water
mesh is decide mesh size and refinement
material props you can change sigmas and D's
solver just solves with the inputs from the other 2

Updated for November 20 check in:

The project is a 2D two-group diffusion solver  using cell centered finite difference and power iteration

To use the files you need pthon, numpy, scipy, and matplotlib

To run the program use the command python verify.py
to make changes to the run you want to do, edit n_asemblies at the top of verify and assembly_pitch in mesh to 
change the reactor size. to change refinement level edit the refinement array at the top of verify. if you want to 
use different sigmas
 for different materials edit the material dictionary at the top of material props, it only works for 2 groups and 
 only works for one type of material. 

verify.py produces verification_order.png when you run it

mesh has functions for defining the size of the reactor
material_props assigns the material values
solver does the finite volume and power iteration over the mesh
exact_answer gives the analytical solution
verify runs both the exact solution and the numerical method and compares the error and plots them

for the material properties i used L. A. Hageman & J. B. Yasinsky Nucl. Sci. Eng. 38, 8–32, 1969, Fig. 3, apparently 
this is the common benchmark example for 2 group diffusuion, found it on ans website

currently the solver is only a solver for a 2D square reactor with a single material and a vacuum boundary. I used 
dirichlet for the basic bare square comparison to verify the solver is working.  to make it actually work for my 
study theres gonna be some changes but the solver works on its own and is ready for the new inputs. the size of the 
reactor i chose is arbitrary, did not take it from the paper. the solver also still works at n = 64, i did a run but 
it took a minute so i went back to 32.  the changes im gonna need to make include:

1. making the mesh into the circle of sqaures I want
2. adding material dictionaries for water and the vacuum after it
3. making multiple FUEL material dictionaries for the multiple fuel zones
4. making atom fraction calcualtions that feed into the material dictonaries for fuel to actually use enrichment, 
for this one i used the predtermined values from the TWIGL benchmark but it should be easy to just place a function 
at the top of material props that spits values into the material dictionaries with enrichment as an input 
5. figure out how to do the zoning
6. solve for enrichment
7. add the water reflector
8. make the marshark boundary condition, not sure if i can just plug that into the else on line 44 of solver or not
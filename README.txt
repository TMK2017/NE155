README for 2-D diffusion solver
====================================
Credits:
Created by Kelsy Green and Joe Corvino
Final Project
Nuclear Engineering 155
University of California, Berkeley
====================================
Purpose:
The neutron diffusion equation describes the averaged behavior of a large amount of neutrons in a non-multiplying medium. The purpose of this code is to write a numerical solver for the two dimensional neutron diffusion equation in MATLAB. There are vacuum boundary conditions on the bottom and left faces and reflecting boundaries on the top and right faces. We will examine the conceptual mathematical principles underlying the code, the algorithms used in the code, how to use the code, and present our code testing techniques.
====================================
How to use:
-Open MATLAB
-Run script "GC_Master.m", input x,y,D,S,sigma,filename. D,S,sigma can be constants, matrices, or function handles. x and y must be vectors. Excel filename must be .xlsx
====================================
Parameter Definitions:
D = diffusion coefficient 
S = source 
sigma  = absorption cross szection



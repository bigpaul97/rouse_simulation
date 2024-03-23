# rouse_simulation

The MATLAB codes for the Rouse polymer simulation described in the paper can be found here:
Main --- Rouse model simulation with uniform friction environment
MainNUF --- Rouse model simulation with non-uniform friction environment
These simulation codes are written for Rouse polymer simulations in 2D. To adapt to other dimensions, please see Instruction below.

The codes are mostly commented and are ready to use given that one puts in the **loop configuration file** and the **save path**
Below are some important instructions for adjusting parameters, boundary conditions, and etc.

Instructions:
-
1. Commonly adjustable Parameters
   --The friction coefficient (zeta) and bead contour length (b) have suggested values for mouse and yeast based on their preferred bead sizes. 
   --The simulation time interval (det) is determined according to the polymer relaxation time.
   --The total simulation timestep (taup_limit) can be changed as like.
   --The boundary condition code (bc) needs to be specified for different boundary conditions.

2. Physical parameters
   --These parameters are usually not changed.
   --One can change the spring constant for the non-nearest-neighbor spring by changing k2.
   
3. Different dimension
   To adapt the simulation to 1-D or 3-D , a few things need to be changed.
   (a) The friction coefficient and pring constant needs to be modified according to Section IIE in the paper.
   (b) In the "Main loop" section of the code, add in (or remove) one dimension for the loop iteration (dmrouse) and saved positions (y1, y2, z1, z2, dimPhyWalk21, dimPhyWalk31, etc.)
   (c) Add in (or remove) one dimension in the "MSD calculation" of the code.

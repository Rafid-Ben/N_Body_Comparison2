# N-Body_Comparison

Updated and condensed version of the N-Body comparison code.
This version is similar to the code in the N-Body_Comparison_Old repository, but hopefully has more features integrated into it instead of added as an afterthought.

Some things this version does not currently have:
* Fragmentation model
* Correct input data/initialization for complex single emitter simulation with different particle types, neutrals, etc.
* Neutral interactions
* Saving of neutral information
* Symmetry codes
* pthread_barrier for multipole method codes (This is implemented in the old version and offers better performance since new threads don't have to be created and destroyed for every function call in the multipole method codes, but it appears there are issues to getting it to work on MAC)

I haven't had much time to test every aspect of the code, so there could be some bugs that need to be fixed and updates that need to be made.
Namely, the writing of data to the timing, L2-norm, energy conservation files hasn't been tested.
Additionally, not all major parameters/conditions for the vectorized code, the BH code, and the multipole method codes have been tested, so it might be worthwhile to run some quick test cases comparing the force calculation method to the direct force method for a new simulation, as the direct force method should hopefully be bug-free and can be used as a ground-truth baseline.

The Barnes-Hut method, setting the domain toggle to DomainToggle::cubicCells, will likely be the easiest version to achieve the fastest speeds.
The fmm codes, especially the vectorized fmm codes, can achieve faster speeds, but it will likely take some work to find optimal parameters for how the cell structure is created for these codes to perform optimally, which will likely differ from computer to computer.

The test_cases.cpp file gives some examples of how to set up and call the code.
The functions are similar to the functions suffixed by \_paper in the test_cases.cpp file in the old version of the code.
I've only tested the compareForceCalculationMethodsAcceleration function so far for a fairly short simulation.

.. module:: source.separation_optimizer

Interface Separation Optimizer
==============================

Class for determining the lowest energy separation at the interface
by finding a minimum of the energy function with respect to the Z-axis.

The solution involves the following steps:

1) The total energy is calculated at a starting guess separation.

2) The separation is changed by a starting step size and the total energy
   is calculated again.

3) If the new energy is lower, then another step is taken in the same
   direction.  If the new energy is greater, then the step size is reduced
   and the direction is reveresed.

4) Step 3 is repeated until:
    
    - the max number of steps is reached
    
    - the step size drops below the provided convergence criteria
    
    - two consecutive steps result in the same energy

The final structure with the lowest energy is returned and the starting
guess is changed to the final separation.

.. tip::
   The method for finding the minimum is designed to find the local 
   minimum near the initial guess and a range dependent on the initial
   step size.  Because of this, is you have trouble getting caught in
   local minima that are not the global minimum, then increase the
   initial step size.

____________

.. autoclass:: SeparationOpt

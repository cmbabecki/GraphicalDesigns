# GraphicalDesigns
Code to compute graphical designs for regular graphs, and a host of example graphs to play around with.  

Regular graphs: The files entitled GDkScript(type).m allow the user to compute k-graphical designs for an input graph running through each valid choice of k of a certain type. The program for combinatorial designs is an integer program that will find the minimum size of a k-graphical design for each k.  The positive and arbitrary weight programs are LP relaxations and hence not guaranteed to find a minimum size design, or even a support minimal design.  The baked-in eigenspace ordering of FindEigenspacesNumeric.m uses the frequency ordering of AD^{-1}, which is introduced in https://arxiv.org/abs/1803.02235 - one can tweak this to use other orderings if desired. Pretty much all common graph operators have the same eigenspaces for regular graphs, with eigenvalues differing by a linear transformation.

To see some worked out examples, go to https://sites.math.washington.edu/~GraphicalDesigns/.  

For a brief overview of the math, see my "What is...?" article in the AMS notices \https://www.ams.org/journals/notices/202209/rnoti-p1571.pdf .

Irregular graphs: The circuits folder contains a program FindCircuits which works for any undirected simple graph. It is pretty straightforward to tweak the scripts for use on positively weighted graphs or with other graph operators, but I haven't done that in a presentable way yet. 

# GraphicalDesigns
Code to compute graphical designs for regular graphs, and a host of example graphs to play around with.  

Regular graphs: The files entitled GDkScript(type).m allow the user to compute k-graphical designs for an input graph running through each valid choice of k of a certain type. The program for combinatorial designs is exact, in that it will find the minimum size of a k-graphical design for each k.  The positive and arbitrary weight programs are relaxations and hence not guaranteed to find a minimum size design.   

To see some worked out examples see https://sites.math.washington.edu/~GraphicalDesigns/.  For a brief overview of the math, see my "What is...?" article in the AMS notices \https://www.ams.org/journals/notices/202209/rnoti-p1571.pdf .

Irregular graphs: The circuits folder contains a program FindCircuits which works for any undirected simple graph. The 

# KnapsackCuts

In a toy example with 7 items, capacity 23, weight 3, 4, 5, 6, 7, 8, and 9, and a relaxation given by x:

n = 7\
b = 23\
a = \[3, 4, 5, 6, 7, 8, 9\]\
x = \[0.0, 1.0, 0.0, 0.5, 1.0, 0.0, 1.0\]

...the function call:

separate_knapsack_cut(n, b, a, x, SepParams)

...with parameters:

SepParams = Dict(\
   "minimum_cut_violation" => 1e-4,\
   "separation_point_precision" => 1e-7,\
   "coefficient_tolerance" => 1e-7\
)

...returns the knapsack cut:

CutData(\[1.0, 1.0, 1.0, 1.0, 1.0\], \[2, 4, 5, 6, 7\], 3.0)

...meaning that:

x\[2\] + x\[4\] + x\[5\] + x\[6\] + x\[7\] <= 3

...is a valid violated knapsack cut.


\ Max-k-cluster (k = 5) problem for the graph:
\
\ 7 12
\ 1 2 1
\ 1 3 1
\ 1 5 1
\ 2 5 1
\ 2 6 1
\ 3 4 1
\ 3 5 1
\ 3 6 1
\ 4 6 1
\ 4 7 1
\ 5 6 1
\ 6 7 1
\
\ Note that we have added all the product constraints
\ to get a tighter semidefinite bound.

max

x1*x2 + x1*x3 + x1*x5 + x2*x5 + x2*x6 + x3*x4 + 
x3*x5 + x3*x6 + x4*x6 + x4*x7 + x5*x6 + x6*x7

st

  x1 + x2 + x3 + x4 + x5 + x6 + x7 = 5
  x1*x2 + x1*x3 + x1*x4 + x1*x5 + x1*x6 + x1*x7 - 4 x1 = 0
  x2*x1 + x2*x3 + x2*x4 + x2*x5 + x2*x6 + x2*x7 - 4 x2 = 0
  x3*x1 + x3*x2 + x3*x4 + x3*x5 + x3*x6 + x3*x7 - 4 x3 = 0
  x4*x1 + x4*x2 + x4*x3 + x4*x5 + x4*x6 + x4*x7 - 4 x4 = 0
  x5*x1 + x5*x2 + x5*x3 + x5*x4 + x5*x6 + x5*x7 - 4 x5 = 0
  x6*x1 + x6*x2 + x6*x3 + x6*x4 + x6*x5 + x6*x7 - 4 x6 = 0
  x7*x1 + x7*x2 + x7*x3 + x7*x4 + x7*x5 + x7*x6 - 4 x7 = 0

bin

 x1 x2 x3 x4 x5 x6 x7

end

\ Max-Cut problem for the graph:
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
\ Note that we have fixed x7 = 0 in the model below.

max

  - 2 x1*x2 - 2 x1*x3 - 2 x1*x5 - 2 x2*x5 - 2 x2*x6 
  - 2 x3*x4 - 2 x3*x5 - 2 x3*x6 - 2 x4*x6 - 2 x5*x6
  + 3 x1 + 3 x2 + 4 x3 + 3 x4 + 4 x5 + 5 x6

bin

 x1 x2 x3 x4 x5 x6

end

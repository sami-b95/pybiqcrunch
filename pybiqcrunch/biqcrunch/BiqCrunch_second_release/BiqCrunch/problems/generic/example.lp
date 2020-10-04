Maximize
  obj:   20 z1 * z3 + 26 z1 * z4 + 23 z2 * z3
        + 8 z2 * z5 + 32 z3 * z4 + 13 z4 * z5

Subject To
  c1: z1 + z2 + z3 + z4 + z5 = 3
  c2: 12 z1 * z3 + 24 z1 * z4 + 14 z2 * z3 + 
      16 z2 * z5 + 28 z3 * z4 + 12 z4 * z5 <= 30

Binary
  z1 z2 z3 
  z4 z5

End

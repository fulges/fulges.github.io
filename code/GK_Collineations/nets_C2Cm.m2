restart

A = QQ[a_0..a_2]
B = QQ[b_0..b_1]
C = QQ[c_0..c_5]

ALL = A**B**C

aa = sub(vars(A),ALL)
bb = sub(vars(B),ALL)
cc = sub(vars(C),ALL)

M_13 = (matrix{{b_0,b_1,0},{0,0,b_0},{0,0,b_1}})
M_14 = (matrix{{b_0,0,0},{0,b_0,0},{0,0,b_1}})
M_15 = (matrix{{b_0,b_1,0},{0,b_0,0},{0,0,b_0}})
M_16 = (matrix{{b_0,b_1,0},{0,b_0,b_1},{0,0,b_0}})
M_17 = (matrix{{b_0,b_1,0},{0,b_0,0},{0,0,b_1}})
M_18 = (matrix{{b_0,0,0},{0,b_0+b_1,0},{0,0,b_0}})
M_19 = (matrix{{b_0,b_1,0,0},{0,b_0,b_1,0},{0,0,0,b_0}})
M_20 = (matrix{{b_0,b_1,0,0},{0,0,b_0,0},{0,0,0,b_0}})
M_21 = (matrix{{b_0,b_1,0,0},{0,0,b_0,b_1},{0,0,0,b_0}})
M_22 = (matrix{{b_0,b_1,0,0},{0,0,b_0,0},{0,0,0,b_1}})
M_23 = (matrix{{b_0,b_1,0,0},{0,b_0,b_1,0},{0,0,b_0,b_1}})
M_24 = (matrix{{b_0,b_1,0,0,0},{0,0,b_0,b_1,0},{0,0,0,0,b_0}})
M_25 = (matrix{{b_0,b_1,0,0,0},{0,b_0,b_1,0,0},{0,0,0,b_0,b_1}})
M_26 = (matrix{{b_0,b_1,0,0,0,0},{0,0,b_0,b_1,0,0},{0,0,0,0,b_0,b_1}})


for j from 13 to 26 list (
    TA = diff(transpose bb,aa*M_j);
    baselocus = minors(2,TA);
    muT = mingens baselocus;
    db = codim baselocus;
    degb = if db == 2 then degree baselocus else "NaN";
    degradb = if db == 2 then degree radical baselocus else "NaN";
    (j,numcols muT, db,degb, degradb))    
-- the table records:
-- the index j
-- the dimension of the system of quadrics
-- the codimension of the base locus
-- if the codimension of the base locus is 2, it records its degree and the degree of its radical
netList oo



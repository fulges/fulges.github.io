restart

-- defining the 6-th cyclotomic extension 
F = QQ[t]/(t^2-t+1)

-- define polynomial ring of variables x_(i,j) 
R = F[x_(0,0)..x_(2,2)]

-- a is a primitive 6th root of 1 and b is its inverse
a = t
b = 1-t

-- define the linear forms of the Waring decomposition
L_1 = -a*x_(0,0) -1/3 *x_(1,1) + b*x_(2,2)
L_2 = -b*x_(0,0) -1/3 *x_(1,1) + a*x_(2,2)
L_3 = -b*x_(0,0) +1/3 *b*x_(1,1) + b*x_(2,2)

L_4 = -x_(0,0) - b*x_(1,2) - 1/3 *a*x_(2,1) 
L_5 = b*x_(0,0) + x_(1,2) - 1/3 *a*x_(2,1) 
L_6 = a*x_(0,0) - a*x_(1,2) - 1/3 *a*x_(2,1) 

L_7 = 1/3 * b* x_(0,1) - a*x_(1,0) + x_(2,2) 
L_8 = 1/3 * b* x_(0,1) - b*x_(1,0) -b* x_(2,2) 
L_9 = 1/3 * a* x_(0,1) - b*x_(1,0) + x_(2,2) 

L_10 = -1/3*a*x_(0,1) +b* x_(1,2) -  x_(2,0)
L_11 = -1/3*b*x_(0,1) +a*   x_(1,2) -   x_(2,0)
L_12 = 1/3*x_(0,1)  -  x_(1,2) -  x_(2,0)

L_13 = x_(0,2) - x_(1,0) - 1/3* x_(2,1)
L_14 = x_(0,2) + b* x_(1,0) + 1/3 *a* x_(2,1)
L_15 = x_(0,2)  + a* x_(1,0) + 1/3*b* x_(2,1)

L_16 = b*x_(0,2) - 1/3 *a* x_(1,1) + x_(2,0) 
L_17 = b*x_(0,2) - 1/3 *b* x_(1,1) -b* x_(2,0) 
L_18 = a*x_(0,2) - 1/3 *b* x_(1,1) + x_(2,0) 


-- define the 3x3 matrix whose entries are x_(i,j)
X = transpose genericMatrix(R,3,3)

-- check that indeed the linear forms above give a decomposition of the determinant.
-- the answer to this command is 0
sum(1..18,j-> L_j^3) - 6*det(X)


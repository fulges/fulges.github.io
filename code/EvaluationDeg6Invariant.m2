-- EvaluationDeg6Invariant.m2

-- the script evaluates the degree 6 invariant described in Luque-Thibon 2003
-- on a point of the tensor network variety TNS(C4,(2,2,2,2),(2,2,2,2))


-- set up the vector spaces as the components of degree one in the polynomial rings
KK = QQ
for i from 1 to 4 do (
    V_i = KK[v_(i,0),v_(i,1)];
    H_i = KK[h_(i,0,(0,0))..h_(i,1,(1,1))])

VVV = V_1**V_2**V_3**V_4; -- V spaces
HHH = H_1**H_2**H_3**H_4; -- Hom spaces


ALL = VVV ** HHH

-- normalization options: the first matrices in the X1 and X3 linear maps are set to be
-- rank one matrices with a single nonzero entry.
normalization = flatten flatten flatten apply(2,k-> (apply(2,i-> apply(2,j-> (
		    if (i,j) == (0,0) then (h_(2*k+1,0,(i,j))=>1) else 
		    h_(2*k+1,0,(i,j)) => 0)))))

-- define the homomorphisms as 2x2 matrices whose entries are linear forms in the v variables
MMM = toList apply(1..4, k-> matrix(apply(2,i-> apply(2,j-> (
		    h_(k,0,(i,j))*v_(k,0) + h_(k,1,(i,j))*v_(k,1))))))

-- the normalized graph tensor
Tmps =sub( trace(product(MMM)),normalization)

-- the 2x2 matrix image of the flattening
mat2x2 = diff(sub(basis({1,0,0,0},VVV),ALL), diff(transpose(sub(basis({0,0,1,0},VVV),ALL)),Tmps))

-- the determinant of the 2x2 matrix, regarded as a bi-quadratic polynomial in the V2 and V4 space
dMat2x2 = det(mat2x2);

-- the 3x3 matrix of the bi-quadratic form
mat3x3 = diff(sub(basis(2,V_2),ALL), transpose(diff(sub(basis(2,V_4),ALL),dMat2x2)));
mat3x3 = sub(mat3x3,HHH);

-- the evaluation of the invariant; this verifies it vanishes identically
invariant = det(mat3x3)


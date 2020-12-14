-- DimensionLowerBound.m2
-- Computes the dimension of the tensor network variety TNS(C4,(m,m,m,m),nn)

-- input the number of vertices d, the bond dimension m and the three local dimensions nn
d = 4
m = 2
nn = (2,2,2,2)

-- define the vector spaces of interest: they are regarded as degree one components of a polynomial ring
WWW = QQ;
VVV = QQ;
ALL = QQ;
for i from 1 to d do (
    V_i = QQ[v_(i,0)..v_(i,nn_(i-1)-1)];
    W_i = QQ[w_(i,0,0)..w_(i,m-1,m-1)];
    VVV = VVV ** V_i;
    WWW = WWW ** W_i;
    ALL = ALL ** W_i **V_i)
use(ALL)


-- define generic matrices in WWW, i.e., the local components of the graph tensor
for i from 1 to d do (
   ww_i = sub(transpose genericMatrix(W_i,m,m),ALL))    

-- the graph tensor
T = trace product(1..d, i-> ww_i)

-- a random point in Hom(W1,V1) + Hom(W2,V2) + Hom(W3,V3)
-- the rank of the differential of the parametrization map at randHom 
-- will provide a lower bound on dim TNS

randHom = flatten flatten toList apply(1..d, k->(
	apply(m,i-> apply(m,j-> w_(k,i,j) => sub(random(1,V_k),ALL)))));

-- compute the image of the differential
-- LL will be a list of elements of multidegree (1,1,1), 
-- which are to be interpreted as elements of V1 \otimes V2 \otimes V3
-- generating the image of the differential of the parametrization map

LL = flatten for i from 1 to d list (
    print i;
    wi = sub(vars(W_i),ALL) ;
    vv = sub(vars(V_i),ALL) ;
    flatten entries (sub( (vv ** diff(wi, T)),randHom)));

-- we obtain the rank of the differential computing the length of 
-- a minimal set of generators of its image
minGen = mingens (ideal LL);
numcols(minGen)
-- the first argument of the min function in the lower bound from Corollary 1.2
expDim = sum(d, i-> m^2*nn_i) -d+1 - sum(d, i-> m^2-1)  
-- the dimension of the ambient space
ambDim = product(toList nn) 



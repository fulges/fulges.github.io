restart

-- list maximal tight supports (indices 0,1,2)
MaxTightSupps = {{{0, 0, 2}, {0, 1, 1}, {1, 0, 1}, {2, 2, 0}}, 
                 {{0, 0, 2}, {0, 2, 1}, {1, 2, 0}, {2, 1, 1}},
		 {{0, 0, 2}, {0, 1, 1}, {0, 2, 0}, {1, 0, 1}, {2, 1, 0}},
		 {{0, 0, 2}, {0, 1, 1}, {1, 0, 1}, {1, 2, 0}, {2, 1, 0}},
		 {{0, 0, 2}, {0, 1, 1}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}},
		 {{0, 0, 2}, {0, 2, 1}, {1, 1, 1}, {2, 0, 1}, {2, 2, 0}}, 
		 {{0, 0, 2}, {0, 1, 1}, {0, 2, 0}, {1, 0, 1}, {1, 1, 0}, {2, 0, 0}}, 
		 {{0, 0, 2}, {0, 2, 1}, {1, 1, 1}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}}, 
		 {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 1, 1}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}}}
   
-- define polynomial rings: 
-- the space of tensors is the multidegree (1,1,1) of a polynomial ring 
-- in three sets of variables

A = QQ[a_0..a_2]
B = QQ[b_0..b_2]
C = QQ[c_0..c_2]

ABC = A**B**C

-- define the tensors T1...T8 and T9_(-1).
TTT = apply(8, i-> ( MM = MaxTightSupps_i;
	sum(MM, m-> a_(m_0) * b_(m_1) * c_(m_2))))
TTT = TTT| {  sum(MaxTightSupps_8, m -> a_(m_0) * b_(m_1) * c_(m_2)) - 2*a_2*b_1*c_0}


-- compute the dimension of their symmetry lie algebra in (gl_3)^{x3}/ C^2.

apply(9, j-> (
	T = TTT_j;
	glAT = (basis({1,0,0},ABC)) ** contract(basis({1,0,0},ABC),T);
	glBT = (basis({0,1,0},ABC)) ** contract(basis({0,1,0},ABC),T);
	glCT = (basis({0,0,1},ABC)) ** contract(basis({0,0,1},ABC),T);
	F = contract( transpose basis({1,1,1},ABC),glAT | glBT | glCT);
	"dim gT_"|toString(j+1)|" ="| toString(27 - rank F - 2) ))
netList oo	


-- introduce disjoint sets of variables to compute kronecker square
A' = QQ[a'_0..a'_2]
B' = QQ[b'_0..b'_2]
C' = QQ[c'_0..c'_2]

Aa = QQ[a_(0,0)..a_(2,2)]
Bb = QQ[b_(0,0)..b_(2,2)]
Cc = QQ[c_(0,0)..c_(2,2)]

-- make a ring with all variables

ALL = A**B**C**A'**B'**C'**Aa**Bb**Cc

-- options to substitute a_i -> a'_i and similarly on other factors
-- relations to substitute a_i*a_j -> a_(i,j) and similarly on other factors

opts = flatten apply(3, i-> {a_i=>a'_i , b_i => b'_i, c_i => c'_i})

Jsub = ideal flatten flatten apply(3, i-> apply(3,j-> {a_i*a'_j - a_(i,j), 
	                         b_i*b'_j - b_(i,j), 
	                         c_i*c'_j - c_(i,j) }))

-- preparing ring for kronecker squares
AABBCC = Aa**Bb**Cc;

-- define kronecker squares
TTTsq = apply(9, j -> (
	T = TTT_j;
	sub( (sub(T,ALL) * sub(sub(T,ALL),opts)) % Jsub , AABBCC)))

-- computing Lie algebra
apply(9, j-> (
	T = TTTsq_j;
	glAT = (basis({1,0,0},AABBCC)) ** contract(basis({1,0,0},AABBCC),T);
	glBT = (basis({0,1,0},AABBCC)) ** contract(basis({0,1,0},AABBCC),T);
	glCT = (basis({0,0,1},AABBCC)) ** contract(basis({0,0,1},AABBCC),T);
	F = contract( transpose basis({1,1,1},AABBCC),glAT | glBT | glCT);
	"dim gT_"|toString(j+1)|" ="| toString(243 - rank F - 2) ))
netList oo
 
 
 -------
 --------
 -- symmetry algebra for T9mu

 restart
 
 KK = QQ[t];
 
 MM =  {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 1, 1}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
 
A = KK[a_0..a_2]
B = KK[b_0..b_2]
C = KK[c_0..c_2]

ABC = A**B**C

T = sum(MM, m -> a_(m_0) * b_(m_1) * c_(m_2)) + (t-1) *a_2*b_1*c_0
glAT = (basis({1,0,0},ABC)) ** contract(basis({1,0,0},ABC),T);
glBT = (basis({0,1,0},ABC)) ** contract(basis({0,1,0},ABC),T);
glCT = (basis({0,0,1},ABC)) ** contract(basis({0,0,1},ABC),T);
F = contract( transpose basis({1,1,1},ABC),glAT | glBT | glCT);

U = F^{1..4,6..25}_{0..16,18..21,23..25}
factor det U


-----
 -- symmetry algebra for T9mu square
 restart
 
 KK = QQ[t];
 
 MM =  {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 1, 1}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
 
A = KK[a_0..a_2]
B = KK[b_0..b_2]
C = KK[c_0..c_2]

ABC = A**B**C

T = sum(MM, m -> a_(m_0) * b_(m_1) * c_(m_2)) + (t-1) *a_2*b_1*c_0

A' = KK[a'_0..a'_2]
B' = KK[b'_0..b'_2]
C' = KK[c'_0..c'_2]

Aa = KK[a_(0,0)..a_(2,2)]
Bb = KK[b_(0,0)..b_(2,2)]
Cc = KK[c_(0,0)..c_(2,2)]

-- make a ring with all variables

ALL = A**B**C**A'**B'**C'**Aa**Bb**Cc
-- options to substitute a_i -> a'_i and similarly on other factors
-- relations to substitute a_i*a_j -> a_(i,j) and similarly on other factors

opts = flatten apply(3, i-> {a_i=>a'_i , b_i => b'_i, c_i => c'_i})

Jsub = ideal flatten flatten apply(3, i-> apply(3,j-> {a_i*a'_j - a_(i,j), 
	                         b_i*b'_j - b_(i,j), 
	                         c_i*c'_j - c_(i,j) }))

-- preparing ring for kronecker squares
AABBCC = Aa**Bb**Cc;

-- define kronecker squares
Tsq = sub( (sub(T,ALL) * sub(sub(T,ALL),opts)) % Jsub , AABBCC)

glAT = (basis({1,0,0},AABBCC)) ** contract(basis({1,0,0},AABBCC),Tsq);
glBT = (basis({0,1,0},AABBCC)) ** contract(basis({0,1,0},AABBCC),Tsq);
glCT = (basis({0,0,1},AABBCC)) ** contract(basis({0,0,1},AABBCC),Tsq);
F = contract( transpose basis({1,1,1},AABBCC),glAT | glBT | glCT);

U = F_{1..9,11..29,31..241};


UU = random(KK^0,KK^239)
for i from 0 to 728 do (
    print i;
    VV = UU||(U^{i});
    if rank VV == numrows(VV) then (UU = VV))

f = det UU;
factor f


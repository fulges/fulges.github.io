restart

A = QQ[a_0..a_2]
B = QQ[b_0..b_1]
C = QQ[c_0..c_5]

Z = QQ[apply(subsets(6,2), jj -> z_jj)]

ALL = A**B**C**Z

aa = sub(vars(A),ALL)
bb = sub(vars(B),ALL)
cc = sub(vars(C),ALL)
zz = sub(vars(Z),ALL)


-- Orbit 19
use(ALL)

M = matrix{{b_0,b_1,0,0},{0,b_0,b_1,0},{0,0,0,b_0}}
T = diff(transpose bb,(aa * M))

muT = apply(subsets(numcols T,2), ss-> (
	z_ss - det (T_ss)))

Bas2 = minors(2,T)
Bas2sat = saturate sub(Bas2,A) -- the base locus is a single point
mingens Bas2
super basis(2,Bas2sat)

C2 = eliminate(flatten entries aa ,ideal muT)
C2 = ideal delete(z_{0,3},flatten entries mingens C2)

(res C2).dd_2 -- the Hilbert-Burch matrix of a cubic scroll in P4


-- Orbit 20 

M = (matrix{{b_0,b_1,0,0},{0,0,b_0,0},{0,0,0,b_0}})
T = diff(transpose bb,(aa * M))

muT = apply(subsets(numcols T,2), ss-> (
	z_ss - det (T_ss)))

Bas2 = minors(2,T)
Bas2sat = saturate sub(Bas2,A) -- the base locus is a line
mingens Bas2
super basis(2,Bas2sat) 

C2 = eliminate(flatten entries aa ,ideal muT) -- the ideal of a linear space



-- Orbit 21

M = (matrix{{b_0,b_1,0,0},{0,0,b_0,b_1},{0,0,0,b_0}})
T = diff(transpose bb,(aa * M))

muT = apply(subsets(numcols T,2), ss-> (
	z_ss - det (T_ss)))


Bas2 = minors(2,T)
Bas2sat = saturate sub(Bas2,A) -- the base locus is a scheme of length 2
mingens Bas2
super basis(2,Bas2sat)

C2 = eliminate(flatten entries aa ,ideal muT)  -- the ideal of a quadric in P3


-- Orbit 22
M = (matrix{{b_0,b_1,0,0},{0,0,b_0,0},{0,0,0,b_1}})
T = diff(transpose bb,(aa * M))

muT = apply(subsets(numcols T,2), ss-> (
	z_ss - det (T_ss)))

Bas2 = minors(2,T)
Bas2sat = saturate sub(Bas2,A) -- the base locus is a set of two points
mingens Bas2
super basis(2,Bas2sat)

C2 = eliminate(flatten entries aa ,ideal muT)  -- the ideal of a quadric in P3

-- Orbit 23
M = (matrix{{b_0,b_1,0,0},{0,b_0,b_1,0},{0,0,b_0,b_1}})
T = diff(transpose bb,(aa * M))

muT = apply(subsets(numcols T,2), ss-> (
	z_ss - det (T_ss)))

Bas2 = minors(2,T)
Bas2sat = saturate sub(Bas2,A) -- the base locus is empty
mingens Bas2 -- the complete linear series of degree 2
super basis(2,Bas2sat) 


-- Orbit 24
M = (matrix{{b_0,b_1,0,0,0},{0,0,b_0,b_1,0},{0,0,0,0,b_0}})
T = diff(transpose bb,(aa * M))

muT = apply(subsets(numcols T,2), ss-> (
	z_ss - det (T_ss)))

Bas2 = minors(2,T)
Bas2sat = saturate sub(Bas2,A) -- the base locus is a single point
mingens Bas2
super basis(2,Bas2sat)

C2 = eliminate(flatten entries aa ,ideal muT)
C2 = ideal delete(null,apply(flatten entries mingens C2, cc -> (
	if degree(cc) == {0,0,0,2} then cc)) )
(res C2).dd_2 -- the Hilbert-Burch matrix of a cubic scroll in P4


-- Orbit 25
M = (matrix{{b_0,b_1,0,0,0},{0,b_0,b_1,0,0},{0,0,0,b_0,b_1}})
T = diff(transpose bb,(aa * M))

muT = apply(subsets(numcols T,2), ss-> (
	z_ss - det (T_ss)))

Bas2 = minors(2,T)
Bas2sat = saturate sub(Bas2,A) -- the base locus is empty
mingens Bas2 -- the complete linear series of degree 2
super basis(2,Bas2sat)


-- Orbit 26

M = (matrix{{b_0,b_1,0,0,0,0},{0,0,b_0,b_1,0,0},{0,0,0,0,b_0,b_1}})
T = diff(transpose bb,(aa * M))

Tmin = apply(subsets(numcols T,2), ss-> (
	z_ss - det (T_ss)))

Bas2 = minors(2,T)
Bas2sat = saturate sub(Bas2,A) -- the base locus is empty
mingens Bas2 -- the complete linear series of degree 2
super basis(2,Bas2sat)

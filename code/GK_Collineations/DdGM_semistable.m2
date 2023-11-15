restart

A = QQ[a_0..a_2]
B = QQ[b_0..b_2]
C = QQ[c_0..c_2]

L = QQ[lam]

ALL = A**B**C**L

aa = sub(vars(A),ALL)
bb = sub(vars(B),ALL)
cc = sub(vars(C),ALL)

-- III
Tt = a_0*b_0*c_0 + a_1*b_1*c_1 + a_2*b_2*c_2 + 
     lam* (a_0*b_1*c_2 + a_2*b_1*c_0 + a_1*b_2*c_0)
M = diff(bb,diff(transpose cc, Tt))
muT = mingens  minors(2,M)
minors(6,diff(transpose sub(basis(2,A),ALL),muT)) 
-- the 2x2 minors span the full linear system of quadrics if lam is nonzero.


families = {"IV.1(i)", "IV.1(ii)", "IV.1(iii)", "IV.2(i)", "IV.2(iii)", "IV.2(vi)", "IV.4(i)", "IV.4(ii)","IV.4(iii)","IV.5(i)","IV.5(ii)", "IV.7(i)","IV.ps", "V"}

T_0 = a_0*b_0*c_0 + a_1*b_1*c_1 + a_2*b_2*c_2 + 
     a_0*b_1*c_2 + a_0*b_2*c_1 + a_1*b_0*c_2 + a_1*b_2*c_0
T_1 = a_0*b_0*c_0 + a_1*b_1*c_1 + a_2*b_2*c_2 + 
     a_1*b_0*c_2 + a_2*b_0*c_1 + a_0*b_1*c_2 + a_2*b_1*c_0
T_2 = a_0*b_0*c_0 + a_1*b_1*c_1 + a_2*b_2*c_2 + 
     a_2*b_1*c_0 + a_1*b_2*c_0 + a_2*b_0*c_1 + a_0*b_2*c_1

T_3 = a_0*b_0*c_0 + a_1*b_1*c_1 + a_2*b_2*c_2 + 
     a_0*b_1*c_2 + a_0*b_2*c_1 + a_1*b_0*c_2
T_4 = a_0*b_0*c_0 + a_1*b_1*c_1 + a_2*b_2*c_2 + 
     a_1*b_0*c_2 + a_2*b_0*c_1 + a_0*b_1*c_2
T_5 = a_0*b_0*c_0 + a_1*b_1*c_1 + a_2*b_2*c_2 + 
     a_2*b_1*c_0 + a_1*b_2*c_0 + a_2*b_0*c_1

T_6 = a_0*b_0*c_0 + a_1*b_1*c_1 + a_2*b_2*c_2 + 
     a_0*b_1*c_2 + a_0*b_2*c_1
T_7 = a_0*b_0*c_0 + a_1*b_1*c_1 + a_2*b_2*c_2 + 
     a_1*b_0*c_2 + a_2*b_0*c_1
T_8 = a_0*b_0*c_0 + a_1*b_1*c_1 + a_2*b_2*c_2 + 
     a_2*b_1*c_0 + a_1*b_2*c_0

T_9 = a_0*b_0*c_0 + a_1*b_1*c_1 + a_2*b_2*c_2 + 
     a_0*b_1*c_2 + a_1*b_2*c_0
T_10 = a_0*b_0*c_0 + a_1*b_1*c_1 + a_2*b_2*c_2 + 
     a_1*b_0*c_2 + a_2*b_1*c_0

T_11 = a_0*b_0*c_0 + a_1*b_1*c_1 + a_2*b_2*c_2 + 
     a_0*b_1*c_2

T_12 = a_0*b_0*c_0 + a_1*b_1*c_1 + a_2*b_2*c_2 

T_13 = a_0*b_1*c_2 + a_1*b_2*c_0 + a_2*b_0*c_1 -  
       (a_1*b_0*c_2 + a_2*b_1*c_0 + a_0*b_2*c_1)

for j from 0 to 13 list (
    TA = diff(bb,diff(transpose cc, T_j));
    baselocus = minors(2,TA);
    muT = mingens baselocus;
    db = codim baselocus;
    degb = if db == 2 then degree baselocus else "NaN";
    degradb = if db == 2 then degree radical baselocus else "NaN";
    {families_j,numcols muT, db,degb, degradb})    
-- the table records:
-- the orbit index in the tables of [DdGM23]
-- the dimension of the system of quadrics
-- the codimension of the base locus
-- if the codimension of the base locus is 2, it records its degree and the degree of its radical
netList oo







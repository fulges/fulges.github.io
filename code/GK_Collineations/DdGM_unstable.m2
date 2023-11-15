restart

A = QQ[a_0..a_2]
B = QQ[b_0..b_2]
C = QQ[c_0..c_2]

ALL = A**B**C

aa = sub(vars(A),ALL)
bb = sub(vars(B),ALL)
cc = sub(vars(C),ALL)

families = {"I.1(i)","I.1(ii)","I.2(i)","I.2(iii)",
            "I.3(i)","I.3(iii)","I.4(i)","I.4(iii)",
	    "I.5(i)","I.5(iii)","I.6(i)","I.6(iii)",
	    "I.7(i)","I.7(ii)","I.8","I.10",
	    "I.11(i)","I.11(ii)","I.14(i)","I.14(iii)"}
	     
T_1 = a_0*b_1*c_2 + a_0*b_2*c_1 + a_1*b_0*c_2 + 
     a_1*b_1*c_1 + a_1*b_2*c_0 + a_2*b_0*c_0
T_2 = a_1*b_0*c_2 + a_2*b_0*c_1 + a_0*b_1*c_2 + 
     a_1*b_1*c_1 + a_2*b_1*c_0 + a_0*b_2*c_0
     
T_3 = a_0*b_1*c_2 + a_0*b_2*c_1 + a_1*b_0*c_2 + 
     a_1*b_1*c_0 + a_1*b_1*c_1 + a_2*b_0*c_0
T_4 = a_2*b_0*c_1 + a_1*b_0*c_2 + a_2*b_1*c_0 + 
     a_0*b_1*c_1 + a_1*b_1*c_1 + a_0*b_2*c_0

T_5 = a_0*b_0*c_2+a_0*b_1*c_1+a_0*b_2*c_0 + 
     a_1*b_0*c_1+a_1*b_1*c_2+a_2*b_0*c_0 
T_6 = a_2*b_0*c_0+a_1*b_0*c_1+a_0*b_0*c_2 + 
     a_1*b_1*c_0+a_2*b_1*c_1+a_0*b_2*c_0

T_7 = a_0*b_0*c_2 + a_0*b_1*c_1 + a_1*b_0*c_1 + 
     a_1*b_1*c_0 + a_2*b_2*c_0 
T_8 = a_2*b_0*c_0 + a_1*b_0*c_1 + a_1*b_1*c_0 + 
     a_0*b_1*c_1 + a_0*b_2*c_2

T_9 = a_0*b_0*c_2 + a_0*b_2*c_0 + a_0*b_2*c_1 + 
     a_1*b_1*c_0 + a_2*b_0*c_1
T_10= a_2*b_0*c_0 + a_0*b_0*c_2 + a_1*b_0*c_2 + 
     a_0*b_1*c_1 + a_1*b_2*c_0

T_11= a_0*b_0*c_2 + a_0*b_1*c_1 + a_1*b_0*c_1 + 
     a_1*b_2*c_0 + a_2*b_1*c_0 
T_12= a_2*b_0*c_0 + a_1*b_0*c_1 + a_1*b_1*c_0 + 
     a_0*b_1*c_2 + a_0*b_2*c_1

T_13= a_0*b_0*c_2 + a_0*b_1*c_1 + a_0*b_2*c_0 + 
     a_1*b_0*c_1 +  a_2*b_1*c_0 
T_14= a_0*b_0*c_2 + a_1*b_0*c_1 + a_2*b_0*c_0 + 
     a_0*b_1*c_1 +  a_1*b_2*c_0 

T_15= a_0*b_0*c_2 + a_0*b_2*c_0 + a_1*b_1*c_1 + 
     a_2*b_0*c_0
     
T_16= a_0*b_0*c_2 + a_0*b_1*c_1 + a_0*b_2*c_0 + 
     a_1*b_0*c_1 + a_1*b_0*c_1 + a_1*b_0*c_1  

T_17= a_0*b_0*c_2 + a_0*b_2*c_0 + a_1*b_0*c_1 +
    a_2*b_1*c_0
T_18= a_0*b_0*c_2 + a_2*b_0*c_0 + a_0*b_1*c_1 +
    a_1*b_2*c_0
    
T_19 = a_0*b_0*c_2 + a_0*b_1*c_0 + a_0*b_2*c_1 +
    a_1*b_0*c_0 + a_2*b_0*c_1
T_20 = a_2*b_0*c_0 + a_0*b_0*c_1 + a_1*b_0*c_2 +
    a_0*b_1*c_0 + a_1*b_2*c_0


for j from 1 to 20 list (
    TA = diff(bb,diff(transpose cc, T_j));
    baselocus = minors(2,TA);
    muT = mingens baselocus;
    db = codim baselocus;
    degb = if db == 2 then degree baselocus else "NaN";
    degradb = if db == 2 then degree radical baselocus else "NaN";
    {families_(j-1),numcols muT, db,degb, degradb})    
-- the table records:
-- the orbit index in the tables of [DdGM23]
-- the dimension of the system of quadrics
-- the codimension of the base locus
-- if the codimension of the base locus is 2, it records its degree and the degree of its radical
netList oo

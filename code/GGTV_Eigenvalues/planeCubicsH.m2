restart

KK = QQ

-- introduce rings:
-- R : the variables of the cubic
-- A : the coefficients of the cubic
-- S : a free variable s on a pencil

R = KK[x_0..x_2]
A = KK[apply(compositions(3,3), ss -> a_ss)]
S = KK[s]

AR = A**R**S
-- initialize vectors for the x and a variables
-- and a general polynomial fgen with variable coefficients
xx = sub(vars(R),AR)
aa = sub(vars(A),AR)

fgen = 6*sum(compositions(3,3), iii-> (
	a_iii * product(3, j-> 1/(iii_j)!*x_j^(iii_j))))

-- the file cubicInvs.m2 contains a hardcoded expression 
-- of the discriminant of plane cubics; we use it to define
-- a handy function that evaluates it.
load "cubicInvs.m2"

Disc = ff -> (
    apply(compositions(3,3), cc -> (
	    mm = product(3, j-> x_j^(cc_j));
	    ee_cc = diff(mm,ff)));
    disc())


-- in each case of the classification of plane cubics
-- we consider a pencil through the normalized element and a generic element.
-- In particular:
-- f is the normalized singular cubic
-- g is the generic cubic, normalized to be in the tangent cone to Disc at f
-- ch is the pencil and C the characteristic polynomial 
-- J is the ideal generated by the coefficients of the characteristic polynomial
-- discg is the discriminant of g

-- in each case, we show that some power of discg belongs to J
-- in particular: requiring the pencil ch to be in mathcal(H)_Disc
-- forces g to be singular: a contradiction

-- in some cases, we further specialize g after inspecting the low degree
-- elements of J. 

-- nodal cubic
f = x_0*x_1^2 - x_2^2*(x_2+x_0)
g = sub(fgen , a_{3,0,0} => 0)
ch = f + s*g

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A);
JJ = mingens J;

discg = sub(Disc(g),A);
discg % J

-- cuspidal cubic
f = x_0*x_1^2 - x_2^3
g = sub(fgen , a_{3,0,0} => 0)
ch = f + s*g

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A);
JJ = mingens J;

JJ_(0,0)

gSpec = sub(g,a_{2,0,1}=>0)
discg = sub(Disc(gSpec),A);
gbTrace = 3
discg^2 % J


-- conic and secant
f = x_0*(x_0^2 + x_1*x_2)
g = sub(fgen , a_{0,3,0} => 0)
ch = f + s*g

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A);
JJ = mingens J;

discg = sub(Disc(g),A);
gbTrace = 3
discg^2 % J


-- conic and tangent
f = x_0*(x_0*x_1 + x_2^2)
g = sub(fgen , a_{0,3,0} => 0)
ch = f + s*g

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A);
JJ = mingens J;
JJ_(0,0)

gSpec = sub(g, a_{0,2,1}=>0)
discg = sub(Disc(gSpec),A);
gbTrace = 3
discg^3 % (J + sub(ideal(a_{0,2,1}),A))


-- triangle 
f = x_0*x_1*x_2
g = sub(fgen , a_{0,3,0} => 0)
ch = f + s*g

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A);
JJ = mingens J;

discg = sub(Disc(g),A);
gbTrace = 3
discg % J


-- asterisk
f = x_0*x_1*(x_0+x_1)
g = sub(fgen , a_{0,0,3} => 0)
ch = f + s*g

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A);
JJ = mingens J;

factor(JJ_(0,0))
gSpec = sub(g, a_{0,1,2}=>0)
discg = sub(Disc(gSpec),A);
gbTrace = 3
discg^2 % (J + sub(ideal(a_{0,1,2}),A))


-- double line and a line
f = x_0^2*x_1
-- component 1: g in H_(0,0,1)
g = sub(fgen , a_{0,0,3} => 0)
ch = f + s*g

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A);
JJ = mingens J;

factor(JJ_(0,0))
Jspec = J + radical ideal (JJ_(0,0));
JJ = mingens Jspec;


discg = sub(Disc(g),A);
gbTrace = 3
discg^3 % (Jspec)


-- component 2: bft in Disc'

---- subcases g(x_0=0) = x_1^3
g = fgen - sub(fgen, x_0=>0) + a_{0,0,3}*x_1^3
ch = f + s*g

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A)

discg = sub(Disc(g),A);
gbTrace = 3
discg^2 % (J)

---- subcases g(x_0=0) = x_2^3
g = fgen - sub(fgen, x_0=>0) + a_{0,0,3}*x_2^3
ch = f + s*g

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A);

discg = sub(Disc(g),A);
gbTrace = 3
discg^2 % (J)

---- subcases g(x_0=0) = x_2 ^2 (x_2+x_1)
g = fgen - sub(fgen, x_0=>0) + a_{0,0,3}*x_2^2*(x_2+x_1)
ch = f + s*g

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A);

discg = sub(Disc(g),A);
gbTrace = 3
discg^2 % (J)

---- subcase bft(x_0=0) = x_2 ^2 x_1
g = fgen - sub(fgen, x_0=>0) + a_{0,0,3}*x_2^2*(x_2+x_1)
ch = f + s*g

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A);

discg = sub(Disc(g),A);
gbTrace = 3
discg^2 % (J)


---- subcase bft(x_0=0) = x_1 ^2 x_2
g = fgen - sub(fgen, x_0=>0) + a_{0,0,3}*x_2^2*(x_2+x_1)
ch = f + s*g

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A);

discg = sub(Disc(bft),A);
gbTrace = 3
discg^2 % (J)



-- triple line
f = x_0^3

-- subcase 1: bft(x_0=0) = x_1^3
g = fgen - sub(fgen, x_0=>0) + a_{0,0,3}*x_1^3
ch = f + s*g

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A);

discg = sub(Disc(g),A);
gbTrace = 3
discg^2 % (J)

-- subcase 2: bft(x_0=0) = x_1^2x_2
g = fgen - sub(fgen, x_0=>0) + a_{0,0,3}*x_1^2*x_2
ch = f + s*g

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A);

discg = sub(Disc(g),A);
gbTrace = 3
discg^3 % (J)

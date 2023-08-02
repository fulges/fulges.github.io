restart

KK = QQ
R = KK[x_0..x_2]
S = KK[s]
A = KK[apply(compositions(3,3), ss -> a_ss)]

ARS = A**R**S

aa = sub(vars(A),ARS)
xx = sub(vars(R),ARS)
Tgen = 6*sum(compositions(3,3), ss -> (
	a_ss * product(3, j-> (1/(ss_j)! * x_j^(ss_j)))))

Disc = f -> (
    apply(compositions(3,3), zz -> (
	    ee_zz = 1/3!*diff(product(3, j-> (x_j^(zz_j))) , f)));
	    load "cubicInvs.m2";
	    return disc)

-- nodal cubic
T = x_0*x_1^2 - x_2^2*(x_2+x_0)
bft = sub(Tgen , a_{3,0,0} => 0)
ch = T + s*bft

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A);
JJ = mingens J;

Tdisc = sub(Disc(bft),A);
Tdisc % J

-- cuspidal cubic
T = x_0*x_1^2 - x_2^3
bft = sub(Tgen , a_{3,0,0} => 0)
ch = T + s*bft

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A);
JJ = mingens J;

JJ_(0,0)

bftSpec = sub(bft,a_{2,0,1}=>0)
Tdisc = sub(Disc(bftSpec),A);
gbTrace = 3
Tdisc^2 % J


-- conic and secant
T = x_0*(x_0^2 + x_1*x_2)
bft = sub(Tgen , a_{0,3,0} => 0)
ch = T + s*bft

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A);
JJ = mingens J;

Tdisc = sub(Disc(bft),A);
gbTrace = 3
Tdisc % J


-- conic and tangent
T = x_0*(x_0*x_1 + x_2^2)
bft = sub(Tgen , a_{0,3,0} => 0)
ch = T + s*bft

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A);
JJ = mingens J;
JJ_(0,0)

bftSpec = sub(bft, a_{0,2,1}=>0)
Tdisc = sub(Disc(bftSpec),A);
gbTrace = 3
Tdisc^3 % (J + sub(ideal(a_{0,2,1}),A))


-- triangle 
T = x_0*x_1*x_2
bft = sub(Tgen , a_{0,3,0} => 0)
ch = T + s*bft

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A);
JJ = mingens J;

Tdisc = sub(Disc(bft),A);
gbTrace = 3
Tdisc % J


-- asterisk
T = x_0*x_1*(x_0+x_1)
bft = sub(Tgen , a_{0,0,3} => 0)
ch = T + s*bft

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A);
JJ = mingens J;

factor(JJ_(0,0))
bftSpec = sub(bft, a_{0,1,2}=>0)
Tdisc = sub(Disc(bftSpec),A);
gbTrace = 3
Tdisc^2 % (J + sub(ideal(a_{0,1,2}),A))


-- double line and a line
T = x_0^2*x_1
-- component 1: bft in H_(0,0,1)
bft = sub(Tgen , a_{0,0,3} => 0)
ch = T + s*bft

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A);
JJ = mingens J;

JJ_(0,0)
Jspec = J + radical ideal (JJ_(0,0));
JJ = mingens Jspec

Tdisc = sub(Disc(bft),A);
gbTrace = 3
Tdisc^3 % (Jspec)


-- component 2: bft in Disc'

---- subcases bft(x_0=0) = x_1^3
bft = Tgen - sub(Tgen, x_0=>0) + a_{0,0,3}*x_1^3
ch = T + s*bft

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A)

Tdisc = sub(Disc(bft),A);
gbTrace = 3
Tdisc^2 % (J)

---- subcases bft(x_0=0) = x_2^3
bft = Tgen - sub(Tgen, x_0=>0) + a_{0,0,3}*x_2^3
ch = T + s*bft

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A)

Tdisc = sub(Disc(bft),A);
gbTrace = 3
Tdisc^2 % (J)

---- subcase bft(x_0=0) = x_2 ^2 (x_2+x_1)
bft = Tgen - sub(Tgen, x_0=>0) + a_{0,0,3}*x_2^2*(x_2+x_1)
ch = T + s*bft

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A)

Tdisc = sub(Disc(bft),A);
gbTrace = 3
Tdisc^2 % (J)

---- subcase bft(x_0=0) = x_2 ^2 x_1
bft = Tgen - sub(Tgen, x_0=>0) + a_{0,0,3}*x_2^2*(x_2+x_1)
ch = T + s*bft

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A)

Tdisc = sub(Disc(bft),A);
gbTrace = 3
Tdisc^2 % (J)



---- subcase bft(x_0=0) = x_1 ^2 x_2
bft = Tgen - sub(Tgen, x_0=>0) + a_{0,0,3}*x_2^2*(x_2+x_1)
ch = T + s*bft

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A)

Tdisc = sub(Disc(bft),A);
gbTrace = 3
Tdisc^2 % (J)



-- triple line
T = x_0^3
-- subcase 1: bft(x_0=0) = x_1^3
bft = Tgen - sub(Tgen, x_0=>0) + a_{0,0,3}*x_1^3
ch = T + s*bft

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A)

Tdisc = sub(Disc(bft),A);
gbTrace = 3
Tdisc^2 % (J)

-- subcase 2: bft(x_0=0) = x_1^2x_2
bft = Tgen - sub(Tgen, x_0=>0) + a_{0,0,3}*x_1^2*x_2
ch = T + s*bft

C = Disc(ch);
J = ideal apply(12, j-> sub(diff(s^j,C),s=>0));
J = sub(J,A)

Tdisc = sub(Disc(bft),A);
gbTrace = 3
Tdisc^2 % (J)


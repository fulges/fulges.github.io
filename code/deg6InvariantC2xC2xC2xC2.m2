restart

KK = QQ
VVV = KK
for i from 1 to 4 do (
    V_i = KK[v_(i,0),v_(i,1)];
    E_i = KK[e_(i,0,(0,0))..e_(i,1,(1,1))])

VVV = V_1**V_2**V_3**V_4;
EEE = E_1**E_2**E_3**E_4

ALL = VVV **EEE

--Tgen = sum((0,0,0,0)..(1,1,1,1),iii-> a_iii * product(1..4, j-> v_(j,iii_(j-1))))


MMM = toList apply(1..4, k-> matrix(apply(2,i-> apply(2,j-> (
		    e_(k,0,(i,j))*v_(k,0) + e_(k,1,(i,j))*v_(k,1))))))

normalization = apply(1..4,k-> matrix(apply(2,i-> apply(2,j-> (
		    if (i,j) == 0 then (e_(k,0,(i,j))=>1
			
opts = flatten flatten flatten toList apply(1..4,p-> apply(2,j-> apply(2,k-> e_(p,0,(j,k)) => 1)))

Tmps = trace(product(MMM));

mat2x2 = diff(sub(basis({1,0,0,0},VVV),ALL), diff(transpose(sub(basis({0,0,1,0},VVV),ALL)),Tmps))

dMat2x2 = det(mat2x2);

mat3x3 = diff(sub(basis(2,V_2),ALL), transpose(diff(sub(basis(2,V_4),ALL),dMat2x2)));

mat3x3 = sub(mat3x3,EEE);

mat3x3_(0,0) * mat3x3_(1,1) * mat3x3_(2,2);

F = sub(mat3x3,opts);

det(F);


 sub(random(1,EEE),ALL)))));


mmm = matrix(MutMat3x3);
Q = mmm * (transpose matrix{{-gg0,gg2,-gg1}});

apply(3,i-> apply(3,j-> (
	    f = gcd(gg1, mat3x3_(i,j));
	    degree f)))

JJJ = ideal(mat3x3^{0})

res oo



---

A1 = matrix apply(2,i-> apply(2,j-> sum(2,k-> e_(1,k,(i,j))*v_(1,k))))
A3 = matrix apply(2,i-> apply(2,j-> sum(2,k-> e_(3,k,(i,j))*v_(3,k))))

rand20 = matrix{{1,0},{0,0}}
rand40 = matrix{{1,0},{0,0}}
rand21 = matrix apply(2,i-> apply(2,j-> (e_(2,1,(i,j)))))
rand41 = matrix apply(2,i-> apply(2,j-> (e_(4,1,(i,j)))))


A10 = A1*rand20
A30 = A3*rand40
A11 = A1*rand21
A31 = A3*rand41

M = det(matrix{{trace(A10*A30), trace(A11*A30)},{trace(A10*A31), trace(A11*A31)}})

mat3x3 = diff(sub(basis(2,V_1),ALL), transpose(diff(sub(basis(2,V_3),ALL),M)));

det(mat3x3);

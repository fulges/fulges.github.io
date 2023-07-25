restart

KK = QQ

S = KK[p_1,p_2,a]
X = KK[x_0..x_2]

SX = S**X 

 -- the linear space TV1 representing the tensor, following Prop. 4.1:
 TV1 = matrix {{x_1 - p_1*x_0 , x_2 - p_2 * x_0 ,0 },
               {x_2 + p_2 * x_0 , p_1*x_1 + (a + p_1^2)*x_0, x_1},
	       { 0 , x_1,-x_0}} 

Mp = mingens minors(2,TV1) 
xx = transpose sub(basis(2,X),SX)

-- the conditions under which the linear series defined by the minors is not full:
depM = radical minors(6,diff(xx, Mp)) 


-- verify that in this case the cubic curve is elliptic curve is singular:
f0 = (det(TV1) ) % depM 
f0guess = ( x_0*x_2^2 - (p_1*x_0 - x_1)^2*(2*p_1* x_0 + x_1)) %depM
f0 - f0guess -- result is 0.

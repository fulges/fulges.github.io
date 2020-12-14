restart

S = {}
delete(null, flatten flatten flatten flatten toList apply(2..4,i1-> toList apply(2..4,i2-> toList apply(2..4, i3-> toList apply(2..4,i4-> (
			iii = {i1,i2,i3,i4};
			s = 1;
			for p from 0 to 3 do (
			    iiir = rotate(p,iii);
			    if member(iiir,S) then s=0);
			for p from 0 to 3 do (
			    iiir = rotate(p,reverse(iii));	      
			    if member(iiir,S) then s=0);
			if s == 1 then S = S|{iii}))))))
	       
netList oo
#S
S = sort(S);
string = ""
for s in S do (
    nn = s;
    load "DimensionLowerBound4Factors.m2";
    ss= toSequence(s);
    string = string | toString (ss) ;
    string = string | " & " ;
    string = string |toString(numcols(minGen) ) | " & ";
    string = string |toString(min{ambDim, expDim})|" \\\\ \n")
string

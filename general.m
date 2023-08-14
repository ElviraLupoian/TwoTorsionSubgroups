// the general programs used in our computations of multi-tangent hyperplanes //

// the following function performs the  Newton-Raphson lifitng   // 
// it's input  is (x,t) - x an approximate vector of lenght t  // 


nr := function(x,t) ;
J0 := Evaluate(J, x) ;
H0 := Evaluate(EE,x) ;
IJ0 := J0^-1 ;
HH0 := Matrix(CC, t, 1, H0) ;
xx := Matrix(CC, t, 1, x) ;
nx := xx - IJ0*HH0 ;
nnx := Eltseq(nx) ;
return nnx ;
end function ;

// the following  4 functions search for short vectors in the appropriate lattice when our initial approximations are real  //
// the first two are built into the third //
// the input of the third function is a complex number  x , approximated  to  p decimal places; and  d is the candidate for the minimal degree of the minimal polynomial//
// the output of the  third function is a (d+1) vector  correspoding the the  coefficients of the polynomial which is a candidate for the minimal polynomial of x //
// the final function transforms the vector  returned by the  third function into a polynomial //



ll := function(x,p,d) ;    
r := [] ;
for i in [1..d-1] do ;
r[i] := [] ;
for j in [1..d+1 ] do ;
if i eq j then 
r[i][j] := 1 ;
else r[i][j] := 0 ;
end if ;
end for ;
end for ;
r1 := [] ;
for i in [1..d+1] do ;
r1[i] := (10^p)*Real(x^(d+1-i)) ;
end for ;
r[d] := [ Round(a) : a in r1 ] ;
r2 := [] ;
for i in [1..d+1] do ;
r2[i] := (10^p)*Real(-CC.1*(x^(d+1-i))) ;
end for ;
r[d+1] := [ Round(a) : a in r2] ; 
return r;
end function ; 

lll := function(x,p,d) ; 
r := ll(x,p,d) ;
MM := Matrix(Integers(), d+1, d+1, r ) ;
M := Transpose(MM) ; 
L := Lattice(M) ;
v := ShortestVector(L) ;
return v ;
end function ;

lllcoord := function(x,p,d) ;
v := Eltseq(lll(x,p,d)) ;
r := ll(x,p,d) ;
r1 := r[d] ;
r2 := r[d+1] ;
a := [] ;
for i in [1..d-1] do ;
a[i] := v[i] ;
end for ;
v1 := &+[ a[i]*r2[i] : i in [1..d-1] ] ;
adC := v[d+1] - v1 ;
ad := adC/r2[d] ;
a[d] := ad; 
v2 := [ a[i]*r1[i] : i in [1..d] ] ;
vv2 := &+v2 ;
ad2C := v[d] - vv2 ;
ad2 := ad2C/r1[d+1] ;
a[d+1] := ad2 ;
return a ;
end function ;

Zu<u> := PolynomialRing(Rationals()) ;

minpoly := function(x,p,d) ;
vv :=  lllcoord(x,p,d) ;
t := &+[vv[i]*u^(d+1-i) : i in [1..d+1] ] ;
return t ;
end function ;


// the following three functions are used to search for relations between the coefficients of the quadritangent planes //
// the first function is built into the second function//
// the input of the second function is  a pair of real approximates pt = [a,b], accuarete to  p decimal  places  and  returns  a vector [c_1, ..c_d+1]  such that : (c_d+1)b + (c_d)a^(d-1)  + .. . + (c_2)a + c_1 = 0 //  
// the  third function converts the  output of the second into an equation //

rell := function(pt,p,d) ;
x := pt[1] ;
y := pt[2] ;
r := [] ;
for i in [1..d-1] do ;
r[i] := [] ;
for j in [1..d+1 ] do ;
if i eq j then 
r[i][j] := 1 ;
else r[i][j] := 0 ;
end if ;
end for ;
end for ;
r1 := [] ;
for i in [1..d] do ;
r1[i] := (10^p)*Real(x^(d-i)) ;
end for ;
r1[d+1] := (10^p)*Real(y) ;
r[d] := [ Round(a) : a in r1 ] ;
r2 := [] ;
for i in [1..d] do ;
r2[i] := (10^p)*Real(-CC.1*(x^(d-i))) ;
end for ;
r2[d+1] := (10^p)*Real(-CC.1*y) ;
r[d+1] := [ Round(a) : a in r2] ; 
return r;
end function ; 

relll := function(x,p,d) ; 
r := rell(x,p,d) ;
MM := Matrix(Integers(), d+1, d+1, r ) ;
M := Transpose(MM) ; 
L := Lattice(M) ;
v := ShortestVector(L) ;
return v ;
end function ;
relllcoord := function(x,p,d) ;
v := Eltseq(relll(x,p,d)) ;
r := rell(x,p,d) ;
r1 := r[d] ;
r2 := r[d+1] ;
a := [] ;
for i in [1..d-1] do ;
a[i] := v[i] ;
end for ;
v1 := &+[ a[i]*r2[i] : i in [1..d-1] ] ;
adC := v[d+1] - v1 ;
ad := adC/r2[d+1] ;
a[d+1] := ad; 
v2 := [ a[i]*r1[i] : i in [1..d-1] ] cat [ a[d+1]*r1[d+1] ]  ;
vv2 := &+v2 ;
ad2C := v[d] - vv2 ;
ad2 := ad2C/r1[d] ;
a[d] := ad2 ;
return a ;
end function ;


re := function(x,p,d);
r := relllcoord(x,p,d) ;
rl := -r[d+1];
v := [ (1/rl)*r[i] : i in [1..d] ] ;
t :=  &+[ v[i]*u^(d-i) : i in [1..d] ] ;
return t;
end function ;


// if the imaginary part of the approximation is small and we think that the algebraic number that we are approximating might be real, we use the following funtions to search for its minimal polynomial 

ll1 := function(x,p,d) ;
r := [] ;
for i in [1..d] do ;
r[i] := [] ;
for j in [1..d+1 ] do ;
if i eq j then 
r[i][j] := 1 ;
else r[i][j] := 0 ;
end if ;
end for ;
end for ;
r1 := [] ;
for i in [1..d+1] do ;
r1[i] := (10^p)*(x^(d+1-i)) ;
end for ;
r[d+1] := [ Round(a) : a in r1 ] ;
return r;
end function ; 

lll1 := function(x,p,d) ; 
r := ll1(x,p,d) ;
MM := Matrix(Integers(), d+1, d+1, r ) ;
M := Transpose(MM) ; 
L := Lattice(M) ;
v := ShortestVector(L) ;
return v ;
end function ;




lll1coord := function(x,p,d) ;
v := Eltseq(lll1(x,p,d)) ;
r := ll1(x,p,d) ;
r2 := r[d+1] ;
a := [] ;
for i in [1..d] do ;
a[i] := v[i] ;
end for ;
v1 := &+[ a[i]*r2[i] : i in [1..d] ] ;
adC := v[d+1] - v1 ;
ad := adC/r2[d+1] ;
a[d+1] := ad; 
return a ;
end function ;

Zu<u> := PolynomialRing(Rationals()) ;

minpoly1 := function(x,p,d) ;
vv :=  lll1coord(x,p,d) ;
t := &+[vv[i]*u^(d+1-i) : i in [1..d+1] ] ;
return t ;
end function ;



rell1 := function(pt,p,d) ;
x := pt[1] ;
y := pt[2] ;
r := [] ;
for i in [1..d] do ;
r[i] := [] ;
for j in [1..d+1 ] do ;
if i eq j then 
r[i][j] := 1 ;
else r[i][j] := 0 ;
end if ;
end for ;
end for ;
r1 := [] ;
for i in [1..d] do ;
r1[i] := (10^p)*Real(x^(d-i)) ;
end for ;
r1[d+1] := (10^p)*Real(y) ;
r[d+1] := [ Round(a) : a in r1 ] ;
return r;
end function ; 

relll1 := function(x,p,d) ; 
r := rell1(x,p,d) ;
MM := Matrix(Integers(), d+1, d+1, r ) ;
M := Transpose(MM) ; 
L := Lattice(M) ;
v := ShortestVector(L) ;
return v ;
end function ;




relll1coord := function(x,p,d) ;
v := Eltseq(relll1(x,p,d)) ;
r := rell1(x,p,d) ;
r2 := r[d+1] ;
a := [] ;
for i in [1..d] do ;
a[i] := v[i] ;
end for ;
v1 := &+[ a[i]*r2[i] : i in [1..d] ] ;
adC := v[d+1] - v1 ;
ad := adC/r2[d+1] ;
a[d+1] := ad; 
return a ;
end function ;


re1 := function(x,p,d);
r := relll1coord(x,p,d) ;
rl := -r[d+1];
v := [ (1/rl)*r[i] : i in [1..d] ] ;
t :=  &+[ v[i]*u^(d-i) : i in [1..d] ] ;
return t;
end function ;


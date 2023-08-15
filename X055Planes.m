// this is the Magma code used to compute the quadritangents to X0(55) stated in X055Quadritangents.txt

// the equations defining the quadritangents that we work with are the following 

Zx<[x]> := FunctionField(Integers(),5) ;
f1 := x[1]*x[3] -x[2]^2 + x[2]*x[4] -x[2]*x[5] -x[3]^2 + 3*x[3]*x[4] + 
x[3]*x[5] -2*x[4]^2 -4*x[5]^2 ;
f2 := x[1]*x[4] -x[2]*x[3] +2*x[2]*x[4] -2*x[2]*x[5] -2*x[3]^2 + 
4*x[3]*x[4] 
+5*x[3]*x[5] -2*x[4]^2 -4*x[4]*x[5] -3*x[5]^2 ;
f3 := x[1]*x[5] -2*x[2]*x[5] -x[3]^2 +2*x[3]*x[4] + x[3]*x[5] -x[4]^2 ;

Za<[a]> := PolynomialRing(Integers(),20) ;
ZX<[X]> := FunctionField(Za,4) ;
F1 := Evaluate(f1, X cat [1]) ;
F2 := Evaluate(f2, X cat [1]) ;
F3 := Evaluate(f3, X cat [1]) ;

G1 := Evaluate(F1, [ X[1], X[2], X[3] , a[1]*X[1] + a[2]*X[2] + a[3]*X[3] 
+ a[4]] ) ;
G2 := Evaluate(F2, [ X[1], X[2], X[3] , a[1]*X[1] + a[2]*X[2] + a[3]*X[3] 
+ 
a[4]] ) ;
G3 := Evaluate(F3, [ X[1], X[2], X[3] , a[1]*X[1] + a[2]*X[2] + a[3]*X[3] 
+ a[4]] ) ;

NG1 := Numerator(G1) ;
NG2 := Numerator(G2) ;
CNG1 := Coefficients(NG1) ;

A := [] ;
for i in [1..10] do ;
A[i] := CNG1[i] ;
end for ;

NG2 := Numerator(G2) ;
CNG2 := Coefficients(NG2) ;

B := [] ;
for i in [1..10] do ;
B[i] := CNG2[i] ;
end for ;

NG3 := Numerator(G3) ;
CNG3 := Coefficients(NG3) ;

C := [] ;
for i in [1..10] do ;
C[i] := CNG3[i] ;
end for ;

g1A := A[8] ;
g2A := A[3]*X[1] + A[6]*X[2] + A[9] ;
g3A := A[1]*X[1]^2 + A[2]*X[1]*X[2] + A[4]*X[1] + A[5]*X[2]^2 + A[7]*X[2] 
+ A[10] ;

g1B := B[8] ;
g2B := B[3]*X[1] + B[6]*X[2] + B[9] ;
g3B := B[1]*X[1]^2 + B[2]*X[1]*X[2] + B[4]*X[1] + B[5]*X[2]^2 + B[7]*X[2] 
+ B[10] ;

g1C := C[8] ;
g2C := C[3]*X[1] + C[6]*X[2] + C[9] ;
g3C := C[1]*X[1]^2 + C[2]*X[1]*X[2] + C[4]*X[1] + C[5]*X[2]^2 + C[7]*X[2] 
+ C[10] ;

alpha := g1A*g2B - g1B*g2A ; 
beta := g1A*g3B - g1B*g3A;

X3 := beta/(-alpha) ;
GG1 := Evaluate(G1, [ X[1], X[2], X3, 1] ) ;
GG2 := Evaluate(G2, [ X[1], X[2], X3, 1] ) ;
GG3 := Evaluate(G3, [ X[1], X[2], X3, 1] ) ;
T1 := Numerator(GG1) ;
T2 := Numerator(GG2) ;
T3 := Numerator(GG3) ;
ZXY<X,Y> := PolynomialRing(Za,2) ;
T1XY := Evaluate(T1, [ X,Y, 1,1] ) ;
T2XY := Evaluate(T2, [ X,Y, 1,1] ) ;
T3XY := Evaluate(T3, [ X,Y, 1,1] ) ;
ZW<W> := PolynomialRing(Za) ;
R := Resultant(T1XY, T3XY, Y ) ;
ZW<W> := PolynomialRing(Za) ;
RW := Evaluate(R, [ W,1]) ;
FRW := Factorization(RW) ;


s2 := FRW[4,1] ;
l := MonomialCoefficient(s2, W^8) ;
h := W^4 + a[5]*W^3 + a[6]*W^2 + a[7]*W + a[8] ;
E := s2 - l*h^2 ;
CE := Coefficients(E) ;
d := Discriminant(h) ;

CE[9] := a[9]*d + 1 ;
CE[10] := a[10]*(a[3] -1) + 1 ;
CE[11] := a[11]*(2*a[3] -1 ) + 1 ;

CE[12] := a[12]*(2*a[1]*a[3] - 2*a[1] - 2*a[3]^2 + 3*a[3] - 2) + 1 ;
CE[13] := a[13]*(2*a[2]*a[3] - 2*a[2] - 2*a[3]^2 + 2*a[3] - 1) + 1 ;
CE[14] := a[14]*(8*a[3]^2 + 2*a[3]*a[4] - 12*a[3] - 2*a[4] + 3) + 1 ;


Za<[a]> := PolynomialRing(Integers(),14) ;
E := [] ;

for i in [1..14] do;
E[i] := Evaluate(CE[i], a cat [ 1 : i in [1..6] ]) ; 
end for ; 

CC := ComplexField(2000) ;
ZZa<[a]> := PolynomialRing(CC,14) ;
EE := [] ;
for i in [1..14] do ;
EE[i] := ZZa ! E[i] ;
end for ;

J := JacobianMatrix(EE) ;
im := CC.1 ;
e := Exp(1) ;
e := CC ! e ;

load "general.m";

// the first orbit of quadritangents is obtained using the following approximation 

a1 := [ -0.131210781858957900000000000000 - 0.335309984329247300000000000000*CC.1, -0.162312458779921700000000000000 + 0.986641789184481000000000000000*CC.1, 0.552993994203229800000000000000 - 0.230274614746367000000000000000*CC.1, 
1.45819980258944230000000000000 + 0.687919421483514100000000000000*CC.1, 4.41284318149945600000000000000 - 20.2156485628203700000000000000*CC.1, -75.4189789156730200000000000000 - 29.9698724095299230000000000000*CC.1, 
-30.2809754987022220000000000000 - 157.398458438756300000000000000*CC.1, 59.5419738873192000000000000000 - 58.2275675302449540000000000000*CC.1, 0.000126151356907789821057632188131 - 6.16904662150548956268235515715*(e^(-5))*CC.1, 
1.76793485626125510000000000000 - 0.910749548424022100000000000000*CC.1, -0.474560942568402600000000000000 - 2.06210797783103200000000000000*CC.1, 1.08713992442503370000000000000 + 0.219888081630770250000000000000*CC.1, 
-0.324948357819111500000000000000 - 1.23272779635115600000000000000*CC.1, 0.367450682643035000000000000000 - 0.0792301721627866500000000000000*CC.1 ];


a1 := [ CC ! b : b in a1 ];
N := [a1] ;
J := JacobianMatrix(EE) ;

// we perform 700 steps of Newton-Raphson 

for i in [2..700] do ;
a := N[i-1] ;
a := [ CC ! b : b in a] ;
N[i] := nr(a,14) ;
end for ;

a := N[700];
CC := ComplexField(2000) ;
p1 := a[1];

// we search for the minimal polynomial of the algebraic number approximated by our initial approximation 

m1 := minpoly(a[1],1700,12) ;
r1 := re([a[1],a[2]], 2000,12) ;
r2 := re([a[1],a[3]], 2000,12);
r3 := re([a[1], a[4]], 2000,12);

// we verify that roots of the above equations are points on our scheme 

SK := SplittingField(m1);
R := Roots(m1, SK) ;
R := [a[1] : a in R ];
r := [ r1, r2, r3] cat [ re([a[1],a[i]], 2000,12) : i in [5..14] ] ;
PT := [ [a] cat [ Evaluate(s,a) : s in r] : a in R ];
EPT := [ [Evaluate(e,a) : e in E ] : a in PT ];
NET := [ a : a in EPT | a ne [ 0 : i in [1..14] ] ] ;
assert NET eq [];

print " One orbit of 4-tangent planes of the form -x[4] + a[1]*x[1] + a[2]*x[2] + a[3]*x[3] + a[4]*x[5] ";
print " the minimal polynomial of a1 is:";
m1;
print " and relations for a[2], a[3], a[4] resp. in terms of a1 are:";
r1;
r2;
r3;

O1 := [m1,r1, r2,r3];
 

// a second orbit of quadritangents is obtained using the following approximation 

a1 := [ 0.0167055338406992700000000000000 + 0.351275065266863900000000000000*CC.1, 0.165999934344472630000000000000 - 0.582145667770257100000000000000*CC.1, 0.677516903369820000000000000000 + 0.209720854867015060000000000000*CC.1, 
0.401048422129924100000000000000 - 0.718341776304055600000000000000*CC.1, 5.60363144749258700000000000000 - 1.45330810861091300000000000000*CC.1, 14.4137482519190510000000000000 - 2.70424592446156400000000000000*CC.1, 
12.8353893551361330000000000000 - 2.42914540949942340000000000000*CC.1, 4.46003348845380300000000000000 - 0.329758521729338300000000000000*CC.1, -0.0117724656383621331106160014145 + 0.0285658292994760582958304963067*CC.1,
2.17926106844116000000000000000 + 1.41724170670566640000000000000*CC.1, -1.17568512200231660000000000000 + 1.38897020035938710000000000000*CC.1, 1.01830260045527530000000000000 - 0.169155235896692200000000000000*CC.1, 
1.67360022436494700000000000000 + 1.46673746335345900000000000000*CC.1, 0.539807026848692900000000000000 + 0.118611163553388550000000000000*CC.1 ];

a1 := [ CC ! b : b in a1 ];
N := [a1] ;
J := JacobianMatrix(EE) ;

// we perform 700 steps of Newton-Raphson 

for i in [2..700] do ;
a := N[i-1] ;
a := [ CC ! b : b in a] ;
N[i] := nr(a,14) ;
end for ;

a := N[700];
CC := ComplexField(2000) ;
p2 := a[1];

// we search for the minimal polynomial of the algebraic number approximated by our initial approximation 

m1 := minpoly(a[1],1700,12) ;
r1 := re([a[1],a[2]], 2000,12) ;
r2 := re([a[1],a[3]], 2000,12);
r3 := re([a[1], a[4]], 2000,12);

// we verify that roots of the above equations are points on our scheme 

SK := SplittingField(m1);
R := Roots(m1, SK) ;
R := [a[1] : a in R ];
r := [ r1, r2, r3] cat [ re([a[1],a[i]], 2000,12) : i in [5..14] ] ;
PT := [ [a] cat [ Evaluate(s,a) : s in r] : a in R ];
EPT := [ [Evaluate(e,a) : e in E ] : a in PT ];
NET := [ a : a in EPT | a ne [ 0 : i in [1..14] ] ] ;
assert NET eq [];

print " One orbit of 4-tangent planes of the form -x[4] + a[1]*x[1] + a[2]*x[2] + a[3]*x[3] + a[4]*x[5] ";
print " the minimal polynomial of a1 is:";
m1;
print " and relations for a[2], a[3], a[4] resp. in terms of a1 are:";
r1;
r2;
r3;

O2 := [ m1, r1, r2, r3];

// the final orbit of quadritangents is obtained using the following approximation 

a1 := [ 0.618049515775443100000000000000 - 0.253264617784191500000000000000*CC.1, 2.01882465080317170000000000000 - 0.218839165616384640000000000000*CC.1, 1.61804951577544300000000000000 - 0.253264617784191500000000000000*CC.1, 
-4.88471543020097600000000000000 + 1.45096848024276380000000000000*CC.1, 5.76893526729169900000000000000 - 8.43881700246932600000000000000*CC.1, -20.7145001990877200000000000000 - 24.9539481579001680000000000000*CC.1, 
-39.9408529608476000000000000000 + 16.5399113467249030000000000000*CC.1, -1.53424683536570700000000000000 + 22.3654841329414540000000000000*CC.1, -0.00515503508351189439291019882114 + 
0.00293662996417205916819183694847*CC.1, -1.38536312330047370000000000000 - 0.567694744449064500000000000000*CC.1, -0.425379888652349700000000000000 - 0.0963585900199587200000000000000*CC.1, 0.603248492333176300000000000000 + 
0.0944232530278480600000000000000*CC.1, 1.85181617853712900000000000000 - 0.609908376044566200000000000000*CC.1, 0.579682511807110200000000000000 + 0.337635857980644900000000000000*CC.1 ];

a1 := [ CC ! b : b in a1 ];
N := [a1] ;
J := JacobianMatrix(EE) ;

// we perform 700 steps of Newton-Raphson 

for i in [2..700] do ;
a := N[i-1] ;
a := [ CC ! b : b in a] ;
N[i] := nr(a,14) ;
end for ;

a := N[700];
p3 := a[1] ;
CC := ComplexField(2000) ;
// we search for the minimal polynomial of the algebraic number approximated by our initial approximation 

m1 := minpoly(a[1],1700,6) ;
r1 := re([a[1],a[2]], 2000,6) ;
r2 := re([a[1],a[3]], 2000,6);
r3 := re([a[1], a[4]], 2000,6);

// we verify that roots of the above equations are points on our scheme 

SK := SplittingField(m1);
R := Roots(m1, SK) ;
R := [a[1] : a in R ];
r := [ r1, r2, r3] cat [ re([a[1],a[i]], 2000,6) : i in [5..14] ] ;
PT := [ [a] cat [ Evaluate(s,a) : s in r] : a in R ];
EPT := [ [Evaluate(e,a) : e in E ] : a in PT ];
NET := [ a : a in EPT | a ne [ 0 : i in [1..14] ] ] ;
assert NET eq [];

print " One orbit of 4-tangent planes of the form -x[4] + a[1]*x[1] + a[2]*x[2] + a[3]*x[3] + a[4]*x[5] ";
print " the minimal polynomial of a1 is:";
print " One orbit of 4-tangent planes of the form -x[4] + a[1]*x[1] + a[2]*x[2] + a[3]*x[3] + a[4]*x[5] ";
print " the minimal polynomial of a1 is:";
m1;
print " and relations for a[2], a[3], a[4] resp. in terms of a1 are:";
r1;
r2;
r3;

O3 := [ m1, r1, r2, r3];


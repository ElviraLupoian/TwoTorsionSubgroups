// this is the Magma code used to compute the quadritangents to X0(42) stated in X042Quadritangents.txt
// the first orbit of quadritangents are represented by points on the following scheme 

Zx<[x]> := FunctionField(Integers(),5);
f1 := x[1]*x[3] - x[2]^2 + x[3]*x[4] ;
f2 := x[1]*x[5] - x[2]*x[5] - x[3]^2 + x[4]*x[5] - x[5]^2;
f3 := x[1]*x[4] - x[2]*x[3] + x[2]*x[4] - x[3]^2 + x[3]*x[4] + x[3]*x[5] - x[4]^2 - 2*x[4]*x[5] ;

Za<[a]> := PolynomialRing(Integers(),16);
ZX<[X]> := FunctionField(Za,4);     
ff1 := Evaluate(f1, X cat [ZX ! 1] ) ;  
ff2 := Evaluate(f2, X cat [ZX ! 1] ) ;
ff3 := Evaluate(f3, X cat [ZX ! 1] ) ; 
F1 := Evaluate(ff1, [ a[1]*X[2] + a[2]*X[3] + a[3]*X[4] + a[4], X[2], X[3], X[4] ] ) ;
F2 := Evaluate(ff2, [ a[1]*X[2] + a[2]*X[3] + a[3]*X[4] + a[4], X[2], X[3], X[4] ] ) ;
F3 := Evaluate(ff3, [ a[1]*X[2] + a[2]*X[3] + a[3]*X[4] + a[4], X[2], X[3], X[4] ] ) ;

f := F2 - (a[1] - 1 )*X[2] ;
f := - f ;
c := ZX ! a[1] - 1 ;
X2 := f/c ;

G1 := Evaluate(F1, [ 1, X2, X[3],X[4] ]);
G2 := Evaluate(F3, [ 1, X2, X[3],X[4] ]);
NG1 := Numerator(G1) ;
DG1 := Denominator(G1) ;

NG2 := Numerator(G2) ;
DG2 := Denominator(G2);

CNG1 := Coefficients(NG1) ;
c1 := CNG1[1] ;
c2 := CNG1[2] ;
c3 := CNG1[3] ;
c4 := CNG1[4] ;
c5 := CNG1[5] ;
c6 := CNG1[6] ;
c7 := CNG1[7] ;
c8 := CNG1[8] ;
c9 := CNG1[9] ;

CNG2 := Coefficients(NG2) ;
b1 := CNG2[1] ;
b2 := CNG2[2] ;
b3 := CNG2[3] ;
b4 := CNG2[4] ;
b5 := CNG2[5] ;
b6 := CNG2[6] ;
b7 := CNG2[7] ;

g1 := c1*X[3]^4 + c2*X[3]^3 + c4*X[3]^2 + c6*X[3] + c9 ;
h1 := c3*X[3]^2 + c5*X[3] + c8 ;

alpha1 := c7 ;
s1 := g1 + h1*X[4] + alpha1*X[4]^2 ;
g2 := b1*X[3]^3 + b3*X[3]^2 + b5*X[3] ;
h2 := b2*X[3]^2 + b4*X[3] + b7 ;

alpha2 := b6 ;
s2 := g2 + h2*X[4] + alpha2*X[4]^2 ;

T1 := -alpha2*h1 +  alpha1*h2 ;
T2 := alpha2*g1 - alpha1*g2 ;

NT1 := Numerator(T1) ;
DT1 := Denominator(T1) ;

NT2 := Numerator(T2) ;  
DT2 := Denominator(T2) ;

X4 := T2/T1 ;
H1 := Evaluate(G1, [ 1,1,X[3], X4] ) ;
H2 := Evaluate(G2, [ 1,1,X[3], X4] ) ;

NH1 := Numerator(H1) ;
DH1 := Denominator(H1) ;
NH2 := Numerator(H2) ;
DH2 := Denominator(H2) ;

ZT<T> := PolynomialRing(Za) ;
H := Evaluate(NH1, [1,1, T, 1] ) ;
l := MonomialCoefficient(H, T^8) ;
h := T^4 + a[5]*T^3 + a[6]*T^2 + a[7]*T + a[8] ;
E := H -l*h^2 ;
eqn := Coefficients(E) ;

d := Discriminant(h) ;
e1 := a[9]*d + 1 ;
e2 := a[10]*(a[1] -1) + 1 ;
e3 := a[11]*(a[3] + 1 ) + 1 ;
e4 := a[12]*(a[1] + a[3]) + 1 ;
eqns := eqn cat [e1, e2, e3,e4] ;

E := eqns ;
s1 := Factorization(E[1]);
E[1] := s1[1,1];
e5 := a[13]*(a[1]*a[3] - 3*a[1] - 3*a[3] + 1 ) + 1 ;
e6 := a[14]*(a[1] + 2*a[2] + a[3] ) + 1  ;
e7 := a[15]*(2*a[1] + a[3] - 1 ) + 1  ;
e8 := a[16]*(a[1]*a[3] - 4*a[1]*a[4] + 5*a[1] - 2*a[3]*a[4] + a[3] + 2*a[4] - 3) + 1 ;


E :=  E cat [e5, e6, e7, e8];
// these are the equations defining the first scheme of quadritangent 

CC := ComplexField(1500) ;
ZZa<[a]> := PolynomialRing(CC,16) ;
EE := [] ;
for i in [1..16] do ;
EE[i] := ZZa ! E[i] ;
end for ;

J := JacobianMatrix(EE) ;
im := CC.1 ;
e := Exp(1) ;
e := CC ! e ;

// the first approximation //
a1 := [-1.336351742793349 + 1.2065776939966408*im,      0.45754607574626277 - 0.7132806181617448*im,    0.8864232983962634 - 0.47463925024866555*im,  >
0.04278761949701573*im, -0.07408466331401076 + 1.9743282067234236*im,   12.521102053351557 - 1.7996587172591652*im,     3.9658619697174182 + 4.5626327>
2.2198720597326522 - 4.212904885173583*im,      -8.508697501490414e-9 - 3.842069688949123*(e^-7)*im,    0.3378980317639447 + 0.1745029314311622*im,   >
0.12543734010788202*im, 0.609519155396034 + 0.9915587858264411*im,      -0.5327134724218096 - 0.15080871600459989*im,   -0.6655863278363146 - 0.993911>
0.24183936394242483 + 0.16825641294065669*im,   -0.14171026332655576 - 0.08845529565327725*im];

a1 := [ CC ! b : b in a1 ];

N := [a1] ;
J := JacobianMatrix(EE) ;
load "general.m";

// we perform 1000 steps of NR //
for i in [2..1000] do ; 
a := N[i-1] ;
a := [ CC ! b : b in a] ;
N[i] := nr(a,16) ;
end for ;

a := N[1000] ;
CC := ComplexField(1500);

// we search for short vectors in an appropriate lattice //

m1 :=  minpoly(a[1], 1500, 8);
r1 := re([a[1],a[2]], 1500, 8) ;
r2 := re([a[1], a[3]], 1500, 8) ;
r3 := re([a[1], a[4]], 1500, 8) ;
print  " One orbit  of 4-tangent planes of the form   -x[1] + a1*x[2] +a2*x[3] + a3*x[4] + a4x[5]" ;


print "The minimal polynomial of a1  is:";
m1;
print "and relations for a2,a3,a4 respectively, in terms of  a1 are:"; 
r1;
r2;
r3;
 

O1  := [ m1, r1, r2, r3];

// the second  scheme of quadritangents that we work with is defined as follows //

Zx<[x]> := FunctionField(Integers(),5);
f1 := x[1]*x[3] - x[2]^2 + x[3]*x[4] ;
f2 := x[1]*x[5] - x[2]*x[5] - x[3]^2 + x[4]*x[5] - x[5]^2;
f3 := x[1]*x[4] - x[2]*x[3] + x[2]*x[4] - x[3]^2 + x[3]*x[4] + x[3]*x[5] - x[4]^2 - 2*x[4]*x[5] ;

Za<[a]> := PolynomialRing(Integers(),16);
ZX<[X]> := FunctionField(Za,4);     
ff1 := Evaluate(f1, X cat [ZX ! 1] ) ;  
ff2 := Evaluate(f2, X cat [ZX ! 1] ) ;
ff3 := Evaluate(f3, X cat [ZX ! 1] ) ; 
F1 := Evaluate(ff1, [ a[1]*X[2] + a[2]*X[3] + a[3]*X[4] + a[4], X[2], X[3], X[4]] ) ;
F2 := Evaluate(ff2, [ a[1]*X[2] + a[2]*X[3] + a[3]*X[4] + a[4], X[2], X[3], X[4]] ) ;
F3 := Evaluate(ff3, [ a[1]*X[2] + a[2]*X[3] + a[3]*X[4] + a[4], X[2], X[3], X[4]] ) ;

f := F2 - (a[1] - 1 )*X[2] ;
f := - f ;
c := ZX ! a[1] - 1 ;
X2 := f/c ;
G1 := Evaluate(F1, [ 1, X2, X[3],X[4] ]);
G2 := Evaluate(F3, [ 1, X2, X[3],X[4] ]);
NG1 := Numerator(G1) ;
DG1 := Denominator(G1) ;
NG2 := Numerator(G2) ;
DG2 := Denominator(G2);

CNG1 := Coefficients(NG1) ;
c1 := CNG1[1] ;
c2 := CNG1[2] ;
c3 := CNG1[3] ;
c4 := CNG1[4] ;
c5 := CNG1[5] ;
c6 := CNG1[6] ;
c7 := CNG1[7] ;
c8 := CNG1[8] ;
c9 := CNG1[9] ;

CNG2 := Coefficients(NG2) ;
b1 := CNG2[1] ;
b2 := CNG2[2] ;
b3 := CNG2[3] ;
b4 := CNG2[4] ;
b5 := CNG2[5] ;
b6 := CNG2[6] ;
b7 := CNG2[7] ;

g1 := c1*X[3]^4 + c2*X[3]^3 + c4*X[3]^2 + c6*X[3] + c9 ;
h1 := c3*X[3]^2 + c5*X[3] + c8 ;

alpha1 := c7 ;
s1 := g1 + h1*X[4] + alpha1*X[4]^2 ;
g2 := b1*X[3]^3 + b3*X[3]^2 + b5*X[3] ;
h2 := b2*X[3]^2 + b4*X[3] + b7 ;

alpha2 := b6 ;
s2 := g2 + h2*X[4] + alpha2*X[4]^2 ;
T1 := -alpha2*h1 +  alpha1*h2 ;
T2 := alpha2*g1 - alpha1*g2 ;NT1 := Numerator(T1) ;
DT1 := Denominator(T1) ;
NT2 := Numerator(T2) ;  
DT2 := Denominator(T2) ;

X4 := T2/T1 ;
H1 := Evaluate(G1, [ 1,1,X[3], X4] ) ;
H2 := Evaluate(G2, [ 1,1,X[3], X4] ) ;

NH1 := Numerator(H1) ;
DH1 := Denominator(H1) ;
NH2 := Numerator(H2) ;
DH2 := Denominator(H2) ;

ZT<T> := PolynomialRing(Za) ;
H := Evaluate(NH1, [1,1, T, 1] ) ;
l := MonomialCoefficient(H, T^8) ;
h := T^4 + a[5]*T^3 + a[6]*T^2 + a[7]*T + a[8] ;
E := H -l*h^2 ;
eqn := Coefficients(E) ;
d := Discriminant(h) ;
e1 := a[9]*d + 1 ;
e2 := a[10]*(a[1] -1) + 1 ;
e3 := a[11]*(a[3] + 1 ) + 1 ;
e4 := a[12]*(a[1] + a[3]) + 1 ;
eqns := eqn cat [e1, e2, e3,e4] ;
E := eqns ;

s1 := Factorization(E[1]);
E[1] := s1[2,1] ;
CNT1 := Coefficients(NT1) ;
t1 := Factorization(CNT1[1])[2,1];
k1 := Factorization(CNT1[2])[2,1];
k2 := Factorization(CNT1[2])[3,1];
t3 := Factorization(CNT1[3])[2,1];
E[13] := a[13]*t1 + 1 ;
E[14] := a[14]*k1 + 1 ;
E[15] := a[15]*k2 + 1 ;
E[16] := a[16]*t3 + 1 ;

CC := ComplexField(1500) ;
ZZa<[a]> := PolynomialRing(CC,16) ;
EE := [] ;
for i in [1..16] do ;
EE[i] := ZZa ! E[i] ;
end for ;


J := JacobianMatrix(EE) ;
im := CC.1 ;
e := Exp(1) ;
// an approxiamtion of a point on this scheme is the following  //
a1 :=   [1.9364169808741407 + 0.4039053131166167*im,    -1.7903167715527235 - 1.625950620227471*im,     -0.15716411214614492 + 0.2600404088677833*im, >
0.3624532321745689*im,  -0.7411908045016656 - 1.5558710219214356*im,    2.501307200160321 - 3.561446905627296*im,       3.2343237448595126 - 0.0816772>
0.21313670860505943 - 0.5055270142029266*im,    9.262496815229144*(e^-5) - 7.176746897424615*(e^-5)*im, -0.9003868624540876 + 0.38836441994689286*im, >
0.3342450285768426*im,  -0.49333729473083304 + 0.1840936676786485*im,   0.1903174923052608 - 0.06221216529039809*im,    0.18118010488839764 - 0.260292>
-0.3189214433991593 + 0.1254057423436788*im,    -0.09964261274141999 + 0.30393542944397648*im] ;
a1 := [ CC ! b : b in a1 ];
N := [a1] ;
J := JacobianMatrix(EE) ;

load "general.m";

// we perform 1000 steps of NR //
for i in [2..1000] do ; 
a := N[i-1] ;
a := [ CC ! b : b in a] ;
N[i] := nr(a,16) ;
end for ;

a := N[1000] ;

CC := ComplexField(1500);
// we search for short vectors in an appropriate lattice //

m1 :=  minpoly(a[1], 1500, 8);
r1 := re([a[1],a[2]], 1500, 8) ;
r2 := re([a[1], a[3]], 1500, 8) ;
r3 := re([a[1], a[4]], 1500, 8) ;

print  " One orbit  of 4-tangent planes of the form   -x[1] + a1*x[2] +a2*x[3] + a3*x[4] + a4x[5]" ;
print "The minimal polynomial of a1  is:";
m1;
print "and relations for a2,a3,a4 respectively, in terms of  a1 are:"; 
r1;
r2;
r3;
 

O2  := [ m1, r1, r2, r3];

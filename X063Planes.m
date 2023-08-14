// this is the magma code used to compute the quadriatanges to X0(63) stated in X0(63)Quadritangents.txt

// the coefficients of the quadritangents are solutions of the following set of equations 

Zx<[x]> := FunctionField(Integers(),5) ;
f1 := x[1]*x[3] - x[2]^2 + x[2]*x[5] - x[3]*x[4] - x[5]^2 ;
f2 := x[1]*x[4] - x[2]*x[3] - x[3]*x[5];
f3 := x[1]*x[5] - x[2]*x[4] - x[3]^2;

Za<[a]> := PolynomialRing(Integers(),20) ;
ZX<[X]> := FunctionField(Za,4) ;
F1 := Evaluate(f1, X cat [1] ) ;
F2 := Evaluate(f2, X cat [1] ) ;
F3 := Evaluate(f3, X cat [1] ) ;

G1 := Evaluate(F1, [ a[1]*X[2] + a[2]*X[3] + a[3]*X[4] + a[4], X[2], X[3], X[4]]) ;
G2 := Evaluate(F2, [ a[1]*X[2] + a[2]*X[3] + a[3]*X[4] + a[4], X[2], X[3], X[4]]) ;
G3 := Evaluate(F3, [ a[1]*X[2] + a[2]*X[3] + a[3]*X[4] + a[4], X[2], X[3], X[4]]) ;

g3 := - X[3]^2 + a[2]*X[3] + a[3]*X[4] + a[4];
g31 := -X[4] + a[1] ;

X2 := g3/(-g31) ;

H1 := Evaluate(G1, [ 1, X2, X[3], X[4]]) ;
H2 := Evaluate(G2, [ 1, X2, X[3], X[4]]) ;
nH1 := Numerator(H1) ;
nH2 := Numerator(H2) ;

ZXY<X,Y> := PolynomialRing(Za,2) ;
nh1 := Evaluate(nH1, [ 1, 1, X, Y]) ;
nh2 := Evaluate(nH2, [ 1, 1, X, Y]) ;

R := Resultant(nh1, nh2, X) ;
FR := Factorization(R)[2,1] ;

ZT<T> := PolynomialRing(Za) ;
H := Evaluate(FR,[ 1,T]) ; 
h := T^4 + a[5]*T^3 + a[6]*T^2 + a[7]*T + a[8] ;
l := MonomialCoefficient(H, T^8) ;
E := H - l*h^2 ; 
CE := Coefficients(E) ;
d := Discriminant(h) ;
CE[9] := a[9]*d + 1 ;
CE[10] := a[10]*a[1] + 1 ;

Za<[a]> := PolynomialRing(Integers(),10) ;
E := [] ;

for i in [1..10] do;
E[i] := Evaluate(CE[i], a cat [ 1 : i in [1..10] ] ) ; 
end for ; 

CC := ComplexField(2000) ;
ZZa<[a]> := PolynomialRing(CC,10) ;
EE := [] ;
for i in [1..10] do ;
EE[i] := ZZa ! E[i] ;
end for ;

J := JacobianMatrix(EE) ;

load "general.m";

im := CC.1 ;
e := Exp(1) ;
e := CC ! e ;

// the first approximation //

a1 := [0.7143533464296229 + 0.35613238129956704*im,     0.8786256402161265 + 0.2651995994997823*im,  0.0284265741394106 - 0.014116479174127233*im,   -0.002722059068912847 - 0.36032223542493763*im, -0.5418350201913729 - 
3.034269878743086*im,   5.84420330315414 + 20.81387303206626*im,        
-7.67468870462932 - 37.89442487731942*im,       -1.1540743391980666 + 
18.38824661801643*im,   -1.2886691073074846*(e^-7) + 
1.721494694208901*(e^-7)*im,    -1.1212032889699048 + 
0.5589625907367356*im];
a1 := [ CC ! b : b in a1 ];
N := [a1] ;
J := JacobianMatrix(EE) ;

// we perform 700 steps of Newton-Raphson 

for i in [2..700] do ;
a := N[i-1] ;
a := [ CC ! b : b in a] ;
N[i] := nr(a,10) ;
end for ;

pt := N[700] ;

// we search for the minimal polynomial of a[1] and relations for the a[2], a[3], a[4] in terms of a[1] 

m1 := minpoly(pt[1],700,8) ;
r1 := re([pt[1],pt[2]], 700, 8) ;
r2 := re([pt[1], pt[3]], 700, 8) ;
r3 := re([pt[1], pt[4]], 700, 8) ;

// we verify that roots of the above equation are points on our scheme 

SK := SplittingField(m1) ;
R := Roots(m1,SK);
R := [ a[1] : a in R ];
r := [ r1,r2, r3] cat [ re([pt[1],pt[i]],700,8) : i in [5..10] ] ;
PT := [ [a] cat [ Evaluate(s,a) : s in r ] : a in R ] ;
EPT := [ [ Evaluate(e,a) : e in E ] : a in PT] ;
NET := [ a : a in EPT | a ne [ 0 : i in [1..10] ]] ;

assert NET eq [] ;

print "One orbit of 4-tangent planes of the form -x[1] + a[1]x[2] + 
a[2]x[3] + a[3]x[4] + a[4]x[5] is as follows:";

print " The minimal polynomial of a[1] is:" ;
m1;

print " and relations for a[2], a[3], a[4] resp. in terms of a[1] are:";
r1;
r2;
r3;

O1 := [ m1, r1, r2, r3];

// We compute a second orbit of planes using the following approxmation  


a1 := [-0.04875698389914024 - 0.7967143359362626*im,    
-0.6689824103483336 + 0.6283123250936408*im,  0.0284265741394106 - 0.014116479174127233*im,   
-0.3106871798919367 + 0.18251849001674914*im,   -2.3568372868337546 + 
1.9863778315173262*im,  -20.947444448490295 - 5.345707990620717*im,     
-7.67468870462932 - 37.89442487731942*im,       16.501725871854553 - 
8.19466561340695*im,    -1.28866910730748469*(e^-7) + 
1.721494694208901*(e^-7)*im,    0.07652584114177508 - 1.2504718264229704*im];


a1 := [ CC ! b : b in a1 ];
N := [a1] ;
J := JacobianMatrix(EE) ;

load "general.m";

for i in [2..700] do ;

a := N[i-1] ;
a := [ CC ! b : b in a] ;
N[i] := nr(a,10) ;
end for ;

pt := N[700] ;


// we search for the minimal polynomial of a[1] and relations for the a[2], a[3], a[4] in terms of a[1]

m1 := minpoly(pt[1],700,8) ;
r1 := re([pt[1],pt[2]], 700, 8) ;
r2 := re([pt[1], pt[3]], 700, 8) ;
r3 := re([pt[1], pt[4]], 700, 8) ;

// we verify that roots of the above equation are points on our scheme

SK := SplittingField(m1) ;
R := Roots(m1,SK);
R := [ a[1] : a in R ];
r := [ r1,r2, r3] cat [ re([pt[1],pt[i]],700,8) : i in [5..10] ] ;
PT := [ [a] cat [ Evaluate(s,a) : s in r ] : a in R ] ;

EPT := [ [ Evaluate(e,a) : e in E ] : a in PT] ;
NET := [ a : a in EPT | a ne [ 0 : i in [1..10] ]] ;

assert NET eq [] ;


print "A second orbit of 4-tangent planes of the form -x[1] + a[1]x[2] +a[2]x[3] + a[3]x[4] + a[4]x[5] is as follows";

print " The minimal polynomial of a[1] is:" ;
m1;

print " and relations for a[2], a[3], a[4] resp. in terms of a[1] are:";
r1;
r2;
r3;

O2 := [ m1, r1, r2, r3];




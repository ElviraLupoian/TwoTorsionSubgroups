// this file test the speed of the implemented lattice reduction 
// the approxiamtion corresponds to the coefficient of the hyperplane, which in the end was not required in our computations 
// the corresponding minimal polynomial has degree 24

// we begin by defining the zero-dimensional scheme 

Zx<[x]> := FunctionField(Integers(),5) ;
f1 := x[1]*x[3] -x[2]^2 + x[2]*x[4] -x[2]*x[5] -x[3]^2 + 3*x[3]*x[4] + x[3]*x[5]-2*x[4]^2 -4*x[5]^2 ;
f2 := x[1]*x[4] -x[2]*x[3] +2*x[2]*x[4] -2*x[2]*x[5] -2*x[3]^2 + 4*x[3]*x[4] +5*x[3]*x[5] -2*x[4]^2 -4*x[4]*x[5] -3*x[5]^2 ;
f3 := x[1]*x[5] -2*x[2]*x[5] -x[3]^2 +2*x[3]*x[4] + x[3]*x[5] -x[4]^2 ;
Za<[a]> := PolynomialRing(Integers(),12) ;
ZX<[X]> := FunctionField(Za,4) ;
F1 := Evaluate(f1, X cat [1]) ;
F2 := Evaluate(f2, X cat [1]) ;
F3 := Evaluate(f3, X cat [1]) ;
G1 := Evaluate(F1, [ X[1], a[1]*X[1] + a[2]*X[3] + a[3]*X[4] + a[4], X[3], X[4]]) ;
G2 := Evaluate(F2, [ X[1], a[1]*X[1] + a[2]*X[3] + a[3]*X[4] + a[4], X[3], X[4]]) ;
G3 := Evaluate(F3, [ X[1], a[1]*X[1] + a[2]*X[3] + a[3]*X[4] + a[4], X[3], X[4]]) ;
NG1 := Numerator(G1) ;
NG2 := Numerator(G2) ;
s1 := - X[3]^2 + 2*X[3]*X[4] + (-2*a[2] + 1)*X[3] - X[4]^2 - 
    2*a[3]*X[4] - 2*a[4];
s2 := ZX ! (2*a[1] -1 ) ;
X1 := s1/s2;
T1 := Evaluate(NG1, [ X1, 1, X[3],X[4]]) ;
T2 := Evaluate(NG2, [ X1, 1, X[3],X[4]]) ;
g1 := Numerator(T1) ;
g2 := Numerator(T2) ;
ZXY<X,Y> := PolynomialRing(Za,2) ;
g1 := Evaluate(g1, [ 1, 1, X,Y]) ;
g2 := Evaluate(g2, [ 1, 1, X,Y]) ;
R := Resultant(g1, g2, Y) ;
ZXYZ<Z> := PolynomialRing(Za) ;
RZ := Evaluate(R, [ Z, 1]) ;
h := Z^4 + a[5]*Z^3 + a[6]*Z^2 + a[7]*Z + a[8] ;
RZF := Factorization(RZ) ;
r := RZF[2,1];
l := MonomialCoefficient(r, Z^8) ;
E := r - l*h^2 ;
CE := Coefficients(E) ;
e9 := a[9]*(Discriminant(h)) + 1;
e10 := a[10]*(2*a[1] - 1 ) + 1;
e11 := a[11]*a[1] + 1 ;
e12 := a[12]*(a[2] + a[3]) + 1 ; 
CE[9] := e9;
CE[10] := e10;
CE[11] := e11;
CE[12] := e12;

Za<[a]> := PolynomialRing(Integers(),12) ;
E := [] ;

for i in [1..12] do;
E[i] := Evaluate(CE[i], a) ; 
end for ; 

CC := ComplexField(2000) ;
ZZa<[a]> := PolynomialRing(CC,12) ;

EE := [] ;
for i in [1..12] do ;
EE[i] := ZZa ! E[i] ;
end for ;

J := JacobianMatrix(EE) ;
im := CC.1 ;
e := Exp(1) ;
e := CC ! e ;

a1 := [0.20370705625860194 + 0.12601289929521822*im,    2.6086646941354554 - 4.678434572156384*im,      2.7388769663106607 + 8.646332076743219*im,      -21.78694937918705 - 
3.200609415489225*im,   -10.75677866940991 - 12.812841080673664*im,     -15.615931279894136 + 75.94632554111494*im,     193.79007420633476 - 169.97577826712057*im,     -101.68868666637773 
+ 473.7118874921864*im, 7.310759366174611*(e^-12) + 4.841150461158379*(e^-12)*im,       1.4290367982737568 + 0.6077669885625143*im,     -3.5503993150476436 + 2.1962720367278306*im,    
-0.12060198290763362 + 0.08948715828938524*im ];
a1 := [ CC ! b : b in a1 ];
N := [a1] ;
J := JacobianMatrix(EE) ;

load "general.m";

for i in [2..700] do ;
a := N[i-1] ;
a := [ CC ! b : b in a] ;
N[i] := nr(a,12) ;
end for ;


pt := N[700] ;



// we now run our code for lattice reduction for  degree 2 to 24 and compare this to the implemented version  
t1 := Cputime() ;
T := [] ;
for i in [1..23] do ;
T[i] := minpoly(pt[1], 700,i+1) ;
 end for ;

m1 :=  T[23] ;
tt1 :=  Cputime(t1) ;

T := [] ;
 t2 := Cputime() ;
for i in [1..23] do ;
T[i] := MinimalPolynomial(pt[1], i+1) ;    
 end for ;                    

tt2 := Cputime(t2) ;

print " To compute the minimal polynomial:";
m1;
print "our lattice reduction takes:";
tt1 ;
print "whilst the implemented version takes:";
tt2 ;


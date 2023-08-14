// This is the MAGMA code used to compute the 2-torsion subgroup of the Jacobian of X072 

Z<U> := PolynomialRing(Integers());
f := U^8 - 12*U^7 + 60*U^6 - 168*U^5 + 348*U^4 - 720*U^3 + 1008*U^2 - 288*U + 144;
K := NumberField(f) ;

// all the quadritangents computed are defined over this field 

OK := Integers(K) ;
PP := PrimesUpTo(50,K) ;
P := PP[3] ;
// this is a prime ideal of norm 5

RR<[x]> := PolynomialRing(Rationals(),5) ;
f1 := x[1]*x[3] -x[2]^2 + x[4]^2 ;
f2 := x[1]*x[4] - x[3]^2 ;
f3 := x[1]*x[5] - x[3]*x[4] -2*x[5]^2 ;
P4 := ProjectiveSpace(RR) ;
X := Curve(P4, [f1,f2,f3] ) ;

// we will reduce the curve modulo P, and use the fact the induced map on the jacobian is an injection when restricted to the torsion subgroup 
F37,pi := ResidueClassField(P) ;
X37 := ChangeRing(X, F37) ;
Cl37, phi, psi := ClassGroup(X37) ;
Z := FreeAbelianGroup(1) ;
degr := hom< Cl37 -> Z | [ Degree(phi(g)) : g in OrderedGenerators(Cl37) ] > ;

J37 := Kernel(degr);   
// these are the points of the Jacobian over the finite field F_5

// we now define the quadritangents computed previously 

R := Roots(f,K) ;
Z<u> := PolynomialRing(K);
S1 := 1/96936*(-67*u^7 + 326*u^6 + 476*u^5 - 3192*u^4 - 3408*u^3 - 
    26712*u^2 + 76104*u + 47184);
S2 := 1/48468*(257*u^7 - 2637*u^6 + 10472*u^5 - 21756*u^4 + 42732*u^3 - 102984*u^2 + 74844*u - 42096);
S3 := 1/6924*(83*u^7 - 800*u^6 + 2924*u^5 - 5760*u^4 + 12696*u^3 - 25608*u^2 + 10512*u + 8928);
coeff1 := [ [-1,R[i,1], Evaluate(S1, R[i,1]), Evaluate(S2, R[i,1]), Evaluate(S3, R[i,1])] : i in [1..8] ] ;

f2 := U^8 - 4*U^7 + 8*U^6 + 8*U^5 - 8*U^4 + 16*U^3 + 32*U^2 - 32*U + 16;
R := Roots(f2,K) ;
Z<u> := PolynomialRing(K);
S1 := 1/32*(-u^7 + 4*u^6 - 6*u^5 - 16*u^4 + 20*u^3 - 40*u + 32);
S2 := 1/32*(-u^7 + 3*u^6 - 2*u^5 - 18*u^4 - 4*u^3 + 20*u^2 - 8*u + 8);
S3 := 1/32*(-u^7 + 4*u^6 - 10*u^5 - 12*u^3 - 32*u^2 + 8*u + 32);
coeff2 := [ [-1,R[i,1], Evaluate(S1, R[i,1]), Evaluate(S2, R[i,1]), Evaluate(S3, R[i,1])] : i in [1..#R] ] ;

f4 := U^8 + 4*U^7 + 20*U^6 + 40*U^5 + 76*U^4 + 80*U^3 + 80*U^2 + 32*U + 16;
R := Roots(f4,K) ;

S1 := 1/184*(15*u^7 + 66*u^6 + 308*u^5 + 668*u^4 + 1076*u^3 + 1152*u^2 + 704*u + 320);
S2 := 1/46*(3*u^7 + 4*u^6 + 34*u^5 - 9*u^4 + 22*u^3 - 64*u^2 - 34*u - 28);
S3 := 1/46*(-5*u^7 - 22*u^6 - 95*u^5 - 192*u^4 - 236*u^3 - 200*u^2 - 112*u + 16);
coeff4 := [ [-1,R[i,1], Evaluate(S1, R[i,1]), Evaluate(S2, R[i,1]), Evaluate(S3, R[i,1])] : i in [1..#R] ] ;


f5 := U^8 + 12*U^6 - 48*U^5 + 156*U^4 - 288*U^3 + 432*U^2 - 288*U + 144 ;
R  := Roots(f5,K) ;
S1 := 1/4884*(-35*u^7 - 2*u^6 - 455*u^5 + 1654*u^4 - 5412*u^3 + 10608*u^2 - 15072*u + 10056);
S2 := 1/2442*(-24*u^7 - 13*u^6 - 312*u^5 + 983*u^4 - 3432*u^3 + 5460*u^2 -7614*u + 1872);
S3 := 1/2442*(-13*u^7 - 24*u^6 - 169*u^5 + 312*u^4 - 1452*u^3 + 312*u^2 - 156*u + 3456);
coeff5 := [ [-1,R[i,1], Evaluate(S1, R[i,1]), Evaluate(S2, R[i,1]), Evaluate(S3, R[i,1])] : i in [1..#R] ] ;


f6 := U^8 - 12*U^7 + 72*U^6 - 24*U^5 + 24*U^4 - 144*U^3 + 288*U^2 - 288*U + 144;
R := Roots(f6,K) ;
S1 := 1/336*(2*u^7 - 25*u^6 + 150*u^5 - 54*u^4 - 288*u^3 - 684*u^2 +  744*u - 456);
S2 := 1/672*(-9*u^7 + 103*u^6 - 594*u^5 - 78*u^4 - 468*u^3 + 1044*u^2 - 1992*u + 984);
S3 := 1/672*(-5*u^7 + 54*u^6 - 294*u^5 - 252*u^4 - 252*u^3 - 72*u^2 - 264*u + 1296);
coeff6 := [ [-1,R[i,1], Evaluate(S1, R[i,1]), Evaluate(S2, R[i,1]), Evaluate(S3, R[i,1])] : i in [1..#R] ] ;

coeff := coeff1 cat coeff2 cat coeff4 cat coeff5 cat coeff6;

// we clear denominatos and any scalars to ensure that our equations are defined over the ring of integers OK

n := #coeff ;
C := [] ;
for i in [1..n] do; 
a1 := coeff[i] ;
a1 := [K ! b  : b in a1 ] ;
da2 := [ Denominator(b) : b in a1 ] ;
lda2 := [ Integers() ! a : a in da2];
L := Lcm(lda2) ;
a3 := [ L*a : a in a1] ;
C[i] := a3 ;
end for ;

CO := [ [ OK ! a : a in b ] : b in C];
Ci := [] ;
for i in [ 1..#CO] do ;
CC := CO[i] ;
C := [ [ Integers() ! a : a in Eltseq(b) ] : b in CC ];
GCi := [ Gcd(a) : a in C ] ;
Ci[i] := GCD(GCi) ;
end for ;

COO := [] ;
for i in [1..#CO] do ;
a := CO[i] ;
ai := [ [ Integers() ! d : d in Eltseq(b)] : b in a ]; 
d  := Ci[i] ;
aid :=  [ [ 1/d*aa : aa in b ] : b in ai] ;
COO[i] := aid ;
end for ;
COK := [ [ OK ! a : a in b ] : b in COO ];
C := COK ;
CO := C ;
ZX<[X]> := PolynomialRing(OK,5) ;
p := [] ;
for i in [1..n] do ;
p[i] := &+[CO[i][j]*X[j] : j in [1..5] ] ;
end for ; 
QT := p ;


// we reduce our equations modulo P 

FX47 := FunctionField(X37) ;
CQT := [ [MonomialCoefficient(a,X[1]),  MonomialCoefficient(a,X[2]) , MonomialCoefficient(a,X[3]) , MonomialCoefficient(a,X[4]), 
MonomialCoefficient(a,X[5])] : a in QT ] ;
F47QT := [ [ pi(s) : s in a ] : a in CQT ] ;

ZU<[U]> := PolynomialRing(F37, 5) ;

F47QT := [ a[1]*U[1] + a[2]*U[2] + a[3]*U[3] + a[4]*U[4] + a[5]*U[5] : a in F47QT];
F47QTT := [ Evaluate(a,[ FX47.1, FX47.2, FX47.3, FX47.4, 1]) : a in F47QT] ;
F47QTT := [ a : a in F47QTT | a ne 0 ] ;

nn := #F47QTT ;
D := [Divisor(F47QTT[i]/F47QTT[1]) : i in [2..nn] ]  ;

DDd := [ Decomposition(a) : a in D ];
DDd1 := [ [ a[i,2] : i in [1..#a] ] :  a in DDd ];
DDd1 := [ [ 1/2*b : b in a ] : a in DDd1 ];
DDd1 := [ [ Integers() ! b : b in a ] : a in DDd1 ];
DDd2 := [ [ a[i,1] : i in [1..#a ] ] : a in DDd ];
divs := [ [ DDd1[i][j]*DDd2[i][j] : j in [1..#DDd1[i] ]  ] : i in [1..#DDd1]];
divs := [ &+a :a in divs ];
H1 := [ psi(a) : a in divs] ;

H := H1 ;
ZN := FreeAbelianGroup(#H) ;
hh := hom< ZN -> J37 | [ a : a in H ] > ;
ihh := Image(hh) ;

print " The 2-torsion subgroup generated by the given orbits of quadritangents is:";

ihh;

assert #ihh eq 2^10 ;
print "Thus we have compute the entire 2-torsion subgroup J0(72)[2]";

// we now take Galois invariants to compute the rational 2-torsion subgroup 

// we found that the Galois group of K has 2 generators, which act on our generators as cpt1 and cpt2 given below 

cpt1 := [ ZN.6 - ZN.7, ZN.3 - ZN.7, ZN.5 - ZN.7, ZN.2 - ZN.7, ZN.4 - ZN.7, - ZN.7, ZN.1 - ZN.7, ZN.13 - ZN.7, ZN.10 - ZN.7, ZN.12 - ZN.7, ZN.15 - ZN.7, ZN.14 - ZN.7, ZN.11 - ZN.7, ZN.9 - ZN.7, ZN.8 - ZN.7, ZN.19 - ZN.7, ZN.22 - ZN.7, ZN.20 - ZN.7, ZN.18 - ZN.7, ZN.16 - ZN.7, ZN.23 - ZN.7, ZN.21 - ZN.7, ZN.17 - ZN.7, ZN.28 - ZN.7, ZN.24 - ZN.7, ZN.31 - ZN.7, ZN.26 - ZN.7, ZN.29 - ZN.7, ZN.25 - ZN.7, ZN.27 - ZN.7, ZN.30 - ZN.7, ZN.34 - ZN.7, ZN.37 - ZN.7, ZN.38 - ZN.7, ZN.33 - ZN.7, ZN.35 - ZN.7, ZN.36 - ZN.7, ZN.39 - ZN.7, ZN.32 - ZN.7 ]; 

conj1 := hom< ZN -> ZN | cpt1>;
mu := hom< ZN -> J37 | [ hh(ZN.i) - hh(conj1(ZN.i)) : i in [1..39]]>;
ker1 := Kernel(mu);
imKer1 := sub<J37 | [hh(k) : k in Generators(ker1)]>;

cpt2 := [ ZN.4 - ZN.3, ZN.7 - ZN.3, -ZN.3, ZN.1 - ZN.3, ZN.6 - ZN.3, ZN.5 - ZN.3, ZN.2 - ZN.3, ZN.14 - ZN.3, ZN.15 - ZN.3, ZN.11 - ZN.3, ZN.10 - ZN.3, ZN.13 - ZN.3, ZN.12 - ZN.3, ZN.8 - ZN.3, ZN.9 - ZN.3, ZN.23 - ZN.3, ZN.20 - ZN.3, ZN.22 - ZN.3, ZN.21 - ZN.3, ZN.17 - ZN.3, ZN.19 - ZN.3, ZN.18 - ZN.3, ZN.16 - ZN.3, ZN.30 - ZN.3, ZN.27 - ZN.3, ZN.29 - ZN.3, ZN.25 - ZN.3, ZN.31 - ZN.3, ZN.26 - ZN.3, ZN.24 - ZN.3, ZN.28 - ZN.3, ZN.35 - ZN.3, ZN.39 - ZN.3, ZN.36 - ZN.3, ZN.32 - ZN.3, ZN.34 - ZN.3, ZN.38 - ZN.3, ZN.37 - ZN.3, ZN.33 - ZN.3];

conj2 := hom< ZN -> ZN | cpt2>;
mu := hom< ZN -> J37 | [ hh(ZN.i) - hh(conj2(ZN.i)) : i in [1..39]]>;
ker2 := Kernel(mu);
imKer2 := sub<J37 | [hh(k) : k in Generators(ker2)]>;

print "The rational 2-torsion subgroup is :";
Rat2tors := imKer1 meet imKer2;
Rat2tors ;

// we now check that this equals the 2-torsion part of the cuspidal subgroup




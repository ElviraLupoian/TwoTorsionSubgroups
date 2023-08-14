// This is the MAGMA code used to compute the 2-torsion subgroup of the Jacobian of the modular curve X0(54) .

// We use field of definition computed in  J054TwoTorsionFOD.m to compute the points of the scheme S and thus the coefficients of the equations defining our tritangent planes.

RR<[x]> := PolynomialRing(Rationals(),4);
 e1 := x[1]^2*x[3] - x[1]*x[3]^2 - x[2]^3 + x[2]^2*x[4] -3*x[2]*x[4]^2 + x[3]^3 + 3*x[4]^3 ;
 e2 := x[1]*x[4] - x[2]*x[3] + x[3]*x[4] ;
 P2 := ProjectiveSpace(RR);
 X := Curve(P2,[e1,e2]);
 ZZ<u> := PolynomialRing(Rationals());
 f1 := u^36 - 24*u^35 + 408*u^34 - 3372*u^33 + 22056*u^32 - 127776*u^31 + 630786*u^30 - 2714424*u^29+ 11261496*u^28 - 41046764*u^27 + 131733144*u^26 - 408414384*u^25 + 1150083423*u^24 - 2814636528*u^23  + 6374836368*u^22 - 13214216088*u^21 + 23592829968*u^20 - 36895248864*u^19+ 51352475964*u^18 - 58328930160*u^17 + 42803579664*u^16 + 3616822728*u^\
15 - 69219558864*u^14 + 126403035264*u^13 - 153084561489*u^12 + 151586318088*u^11 - 131113127592*u^10 + 78540995524*u^9 + 18513274440*u^8 - 121585972992*u^7 + 155356552290*u^6 - 110860906584*u^5 + 57084024120*u^4 - 38472387036*u^3 + 36135508152*u^2 - 21483956688*u + 5029788241;

N := NumberField(f1);  

// we compute the ring of integers ON  of the number field N 
e := EquationOrder(N);
F := Factorization(Discriminant(e));
P := [ a[1] : a in F ];

for p in P do ;
e := pMaximalOrder(e,p);
end for ;
ON := MaximalOrder(e); 

// we reduce our curve modulo P , an ideal in ON of norm 289
//  this reduction map induces an injective map from the K-rational torsion subgroup of J0(54)  to J0(54)(F_289) and use this to verify that we have the entire two-torsion subgroup 

P := 289*ON ;
FP := Factorization(P);
P1 := FP[1,1];
F289,pi := ResidueClassField(P1);
X289 := ChangeRing(X,F289);
Cl289,phi,psi := ClassGroup(X289);
Z := FreeAbelianGroup(1) ;
degr := hom<Cl289 -> Z | [ Degree(phi(g)) : g in OrderedGenerators(Cl289)]>;
 J289 := Kernel(degr) ; // these are the F_289 points of the Jacobian J0(54) 


// we redefine the scheme whose points correspond to tritangent planes and check that the field of definition computed in  J054TwoTorsionFOD.m is in fact the field of definition of the 2-torsion subgroup  
Zxyz<x,y,z>:=PolynomialRing(Integers(),3);
F := x^2*z - x*z^2 - y^3 + y^2 - 3*y +z^3 + 3 ;
G := x - y*z + z ;
Zu<[u]> := PolynomialRing(Integers(),7) ;
ZYZ<Y,Z> := PolynomialRing(Zu,2) ;
f := Evaluate(F, [u[1]*Y+u[2]*Z+u[3] , Y, Z]);
g := Evaluate(G, [u[1]*Y+u[2]*Z+u[3] , Y, Z]);
h := Resultant(f,g,Z);
ZW<W> := PolynomialRing(Zu);
h := Evaluate(h,[W,0]);
l := Coefficient(h,6);
eqns := Coefficients(h - l*(W^3 + u[4]*W^2 + u[5]*W + u[6])^2); 
 p := W^3 + u[4]*W^2 + u[5]*W + u[6] ;
e := Discriminant(p)*u[7] -1 ;
eqns := eqns cat [e] ;
S := Scheme(AffineSpace(Zu),eqns);
SN := BaseChange(S,N);

// the following takes approximately 15 minutes to complete  

A := Points(SN); 
t := #A;

assert t eq 120 ;

print "The field of definition of the 2-torsion subgroup is N ; as expected";

// we clear denominators to ensure that all defining equations are defined over ON and are in minimal form 
A := [ Eltseq(a) : a in  A];
A3 := [ [u[i] : i in [1..3] ] : u in  A];
 A4 :=  [ [-1] cat u : u in  A3];
DA4 := [ [Denominator(aa) : aa in a ] : a in A4 ];
LCMDA4 := [ LCM(a) : a in DA4 ];
CA4 := [ [LCMDA4[i]*a : a in A4[i]] : i in [1..120]]; 
ONCA4 := [ [ ON ! aa : aa in a ] : a in CA4 ];
tf,u := IsPrincipal(P1);
vONCA4 := [ [Valuation(aa,P1) : aa in a ] : a in ONCA4 ];
mVal := [ Minimum(a) : a in vONCA4];
OONCA4 := [ [u^(-mVal[i])*aa : aa in ONCA4[i] ] : i in [1..120] ]; 

// we reduce our coefficients ( and thus the defining equations modulo P 
F289CA4 := [ [pi(aa) : aa in a ] : a in OONCA4 ]; 

// we form our equations defining the tritangent planes over F_289 
TT<x,y,z,w> := PolynomialRing(F289,4) ;
F289TP := [ a[1]*x + a[2]*y + a[3]*z + a[4]*w : a in F289CA4];
NZF289TP := [ a : a in F289TP | a ne 0 ]; 
FX289 := FunctionField(X289);
FFTP := [ Evaluate(a, [ FX289.1, FX289.2, FX289.3,1]) : a in NZF289TP ];

// we now form the divisors ( corresponding to 2-torsion points) corresponding to the tritangent planes 
Db := [ Divisor(FFTP[i]/FFTP[1]) : i in [1..#FFTP]];
DDb := [Decomposition(s) : s in Db];
DDb1 := [ [1/2*a[i,2] : i in [1..#a]] : a in DDb ];
DDb1 := [ [Integers()!s : s in a ] : a in DDb1] ;
DDb2 := [ [a[i,1] : i in [1..#a] ] : a in DDb];
divs := [ [DDb1[i][j]*DDb2[i][j] : j in [1..#DDb1[i]]] : i in [1..#DDb1] ];  
divs := [ &+a : a in divs ];

Z := FreeAbelianGroup(1);
degr := hom<Cl289 -> Z | [ Degree(phi(g)) : g in OrderedGenerators(Cl289)] > ;
JF := Kernel(degr) ;
H := [ psi(a) : a in divs];
H := [JF!s : s in H ];
ZN := FreeAbelianGroup(#H);
h := hom<ZN -> JF | [ a : a in H ]> ;

// the image of h is the subgroup formed by our divisors 

ih := Image(h) ;

assert #ih eq 2^8 ;

print "The 2-torsion subgroup formed using the tritangents  computed is :";
ih;
print "Hence we have computed the entire 2-torsion subgroup of J0(54).";

// we check the fact that the rational 2-torsion subgroup is trivial 

A1 := [ a[1] : a in A];
Aut := Automorphisms(N) ;

T := [] ;
for k in [1..36] do ;
b := Aut[k] ;
cpt := [];
for i in [1..120] do ;
s := [ j : j in [1..120] |  b(A1[i]) eq A1[j] ][1] ;
cpt[i] := s ;
end for ;
cc := [] ;
s := cpt[1] ;
for i in [1..120] do ;
t := cpt[i] ;
cc[i] := ZN.t - ZN.s ;
end for ;
conj := hom< ZN -> ZN | cc >;
mu := hom< ZN -> J289 | [ h(ZN.i) - h(conj(ZN.i)) : i in [1..120]] >;
ker := Kernel(mu) ;
imker := sub< J289 | [ h(k) : k in Generators(ker)]>;
T[k] :=  #imker ;
end for ;


TT := [ i :  i in [1..36] | T[i] eq 1 ];
tt := #TT ;
assert tt ne 0 ;

print "The rational 2-torsion subgroup is trivial.";


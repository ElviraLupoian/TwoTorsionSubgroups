// this is the MAGMA code used to compute the possible field of definion of the  2-torsion subgroup of the Jacobian of the non-hyperelliptic, genus 4 curve X0(54) 
// In the file J054TwoTorsion.m we vefiy that the field of definition computed using this code is in fact the field of definition of the 2-torsion subgroup 
// we construct our 2-torsion points using the tritangent planes to the curve 

// we begin by finding a set of equations whose solutions correspond to coefficients of the defining equations of the tritangents 

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

// note eqns  is the a set of equations defining a zero-dimensional scheme whose points correspond to coefficients of defining equations of tritangent planes to the curve 

S := Scheme(AffineSpace(Zu),eqns);

// we base change this scheme to the finite field F_289, compute its points over this finite field and hensel lift them 

K := QuadraticField(5) ;
OK := Integers(K) ;
P := 17*OK ;
FP := Factorization(P) ;
P := FP[1,1];

F289,phi := ResidueClassField(P);
S2 := BaseChange(S,F289);
A := Points(S2) ;
A := [ [ a@@phi : a in Eltseq(b) ] : b in A ];

Zz<[u]> := PolynomialRing(OK,7);
f1 := Zz ! eqns[1];
f2 := Zz ! eqns[2];
f3 := Zz ! eqns[3];
f4 := Zz ! eqns[4];
f5 := Zz ! eqns[5];
f6 := Zz ! eqns[6];
f7 := Zz ! eqns[7];

H := [f1,f2,f3,f4,f5,f6,f7];
J := JacobianMatrix(H) ;

// we perform the Hensel lifing using the following functions 
// the input of this function is a solution of our equations modulo  p^n and its output is the solutions modolu p^n+1 lifiting the inputed solution 


lift := function(n,A) ;
a := Eltseq(A) ;
Fn,tn := quo<OK | P^(n) > ;
b := [ s@@tn : s in a ];
c := Evaluate(H,b);
d := [ (1/17^n)*s : s in c] ;
e := [ phi(s) : s in d ] ;
Y := Matrix(F289, 7,1, e ) ;
Jc := Evaluate(J,b) ;
s1 := Eltseq(Jc[1]) ;
s2 := Eltseq(Jc[2]) ;
s3 := Eltseq(Jc[3]) ;
s4 := Eltseq(Jc[4]) ;
s5 := Eltseq(Jc[5]) ;
s6 := Eltseq(Jc[6]) ;
s7 := Eltseq(Jc[7]) ;
d1 := [phi(s) : s in s1 ];
d2 := [phi(s) : s in s2 ];
d3 := [phi(s) : s in s3 ];
d4 := [phi(s) : s in s4 ];
d5 := [phi(s) : s in s5 ];
d6 := [phi(s) : s in s6 ];
d7 := [phi(s) : s in s7 ];
s := d1 cat d2 cat d3 cat d4 cat d5 cat d6 cat d7 ;
Jcc := Matrix(F289, 7, 7, s);
Zu<[y]> := PolynomialRing(F289, 7) ;
B := Matrix(Zu, 7,1, [y[1],y[2],y[3],y[4], y[5], y[6], y[7]]) ;
JCc := RMatrixSpace(Zu, 7,7) ! Jcc ;
eqns := Eltseq( (JCc * B ) + Y ) ;
Zzz<[y]> := PolynomialRing(F289,7);
g1 := Zzz ! eqns[1] ;
g2 := Zzz ! eqns[2] ;
g3 := Zzz ! eqns[3] ;
g4 := Zzz ! eqns[4] ;
g5 := Zzz ! eqns[5] ;
g6 := Zzz ! eqns[6] ;
g7 := Zzz ! eqns[7] ;
EQN := [g1,g2,g3,g4 , g5, g6 ,g7 ] ;
S1 := Scheme(AffineSpace(Zzz) , EQN) ;
PT := Points(S1) ;
PPT := [ Eltseq(a) : a in PT ] ;
PTT := [ s@@phi : s in PPT ] ;
NPT := [ [b[1] + (17^n)*s[1], b[2] + (17^n)*s[2], b[3] + (17^n)*s[3], b[4] + (17^n)*s[4], b[5] + (17^n)*s[5] ,b[6] + (17^n)*s[6] ,b[7] + (17^n)*s[7] ] : s in PTT ] ;
return NPT[1] ;
end function ;

// we hensel lift all points of A, to precision 17^200

Approx := [] ;
for i in [1..120] do ;
Approx[i] := [ lift(1, A[i]) ];
end for ;


for j in [1..120] do ;
for i in [2..400] do ;
Approx[j][i] := lift(i, Approx[j][i-1]) ;
end for ;
end for ;


// we now compute candidates for the minimal polynomials of the first coefficient 
// the following function uses lattice reduction to compute such a candidate 
// the input is an approximation a, accurate modulo p^k, and d is a possible degree 

apr := function(x,k,p,d) ; 
ZZ := FreeAbelianGroup(d+1) ;
Z2 := FreeAbelianGroup(2) ;
Z2s := sub< Z2 | [(p^k)*Z2.1,(p^k)*Z2.2] > ;
Q,pi := quo< Z2 | Z2s > ; 
R,r := quo<OK | P^k> ;
a := [ Z2 !  Eltseq((x^i)@@r) : i in [0..d] ] ;
c := [ pi(s) : s in a ] ;
phi := hom< ZZ -> Q | c > ;
K := Kernel(phi) ;
G := { Eltseq(ZZ ! g) : g in Generators(K) } ;
W := StandardLattice(d+1) ;
L := sub<W | [ W ! g : g in G ] > ;
v := ShortestVector(L) ;
II := Index(W,L) ;
Bound := (II)^(1/(d+1));
Bound2 := (1/1000)*Bound ;
Zx<u> := PolynomialRing(Integers());
O := [ v[s+1]*u^s : s in [0..d] ] ;
T := &+[ a : a in O] ;
return T, Length(v),Bound2;
end function ;

// we compute candidates for the minimal polynomials 

M := [] ;

for i in [1..120] do ;
a := Approx[i][400][1] ;
m := apr(a, 400,17,36) ;
M[i] := m;
end for ;


 MM :=[] ;

 for i in [1..120] do ;
 F, pi := quo<OK |P^400>;
 a1 := Approx[i][400][1] ;
 fa := Factorization(M[i]) ;
 s := [ a[1] : a in fa | pi(Evaluate(a[1], a1) ) eq 0 ] ;
MM[i] := s ;
end for ;

t := &*[ #a : a in MM ];

assert t eq 1 ;

MM := [ a[1] : a in MM ] ;

// M is the sequence of candidates for the minimal polynomials of a_1 

MM := Set(MM) ;
MM := SetToSequence(MM) ;


// we verify that all these polynomials factorize completely over the number field defined by the following polynomial 

 
ZZ<u> := PolynomialRing(Rationals());
g  := u^36 - 24*u^35 + 408*u^34 - 3372*u^33 + 22056*u^32 - 127776*u^31 + 630786*u^30 - 2714424*u^29+ 11261496*u^28 - 41046764*u^27 + 131733144*u^26 - 408414384*u^25 + 1150083423*u^24 - 2814636528*u^23  + 6374836368*u^22 - 13214216088*u^21 + 23592829968*u^20 - 36895248864*u^19+ 51352475964*u^18 - 58328930160*u^17 + 42803579664*u^16 + 3616822728*u^\
15 - 69219558864*u^14 + 126403035264*u^13 - 153084561489*u^12 + 151586318088*u^11 - 131113127592*u^10 + 78540995524*u^9 + 18513274440*u^8 - 121585972992*u^7 + 155356552290*u^6 - 110860906584*u^5 + 57084024120*u^4 - 38472387036*u^3 + 36135508152*u^2 - 21483956688*u + 5029788241;

// note that this is the largest degree polynomial in our list of candidates for the minimial polynomials 

K := NumberField(g) ;

// we verify that all polynomials in MM split completely over this number field 

RR := [] ;
t := #MM;

for i in [1..t] do ;
s := #Roots(MM[i], K) ;
d := Degree(MM[i]) ;
n :=  d -s ;
RR[i] := n ;
end for ;

r := [ 0 : i in [1..t] ] ;

assert RR eq r ;

print "A candidate for the field of definition of the 2-torsion subgroup of J0(54) is:";

K;




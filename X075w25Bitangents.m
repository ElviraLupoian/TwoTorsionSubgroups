// This the MAGMA code used to compute the bitangents to the curve X0(75)/w25 ( the modulo curve quotiented by the Atkin-Lehner involution w25).

// We begin by defining a system of equations whose solutions correspond to coefficients of lines which are bitangents to the curve.

Z<x,y,z> := PolynomialRing(Integers(),3);
f := 3*x^3*z - 3*x^2*y^2 + 5*x^2*z^2 -3*x*y^3 -19*x*y^2*z -x*y*z^2 +3*x*z^3 +2*y^4 +7*y^3*z -7*y^2*z^2 -3*y*z^3;
Zz<X,Y> := PolynomialRing(Integers(),2);
F := Evaluate(f,[X,Y,1]);
Zu<[u]> := PolynomialRing(Integers(),4);
ZY<Y> := PolynomialRing(Zu);
h := Evaluate(F,[u[1]*Y + u[2],Y]);
l := MonomialCoefficient(h,Y^4);
H := h - l*(Y^2 + u[3]*Y + u[4])^2 ;
eqns := Coefficients(H);

// eqns is the required system of equations, which define a zero dimensional scheme 

S := Scheme(AffineSpace(Zu),eqns);

// we base change this scheme to the finite field F_289, search for points over this field and then Hensel lift.

K := QuadraticField(7) ; 
OK := Integers(K) ;
P := 17*OK ;
F289,phi := ResidueClassField(P);
S2 := BaseChange(S,F289);
A := Points(S2) ;

A := SetToSequence(A) ;
A := [ Eltseq(a) : a in A ]; 
A := [ [ a@@phi : a in b ] : b in A]; 


Zz<[u]> := PolynomialRing(OK,4);
f1 := Zz ! eqns[1];
f2 := Zz ! eqns[2];
f3 := Zz ! eqns[3];
f4 := Zz ! eqns[4];

H := [f1,f2,f3,f4];
J := JacobianMatrix(H) ;


// the following function will lift a solution of our system modulo p^n to a solution modulo p^n+1 
// its input is a solution A, modulo p^n and n 
// its output is a solution modulo p^n+1 lifting A 


lift := function(n,A) ; 
a := Eltseq(A) ;
Fn,tn := quo<OK | P^(n) > ;
b := [ s@@tn : s in a ];
c := Evaluate(H,b);
d := [ (1/17^n)*s : s in c] ;
e := [ phi(s) : s in d ] ;
Y := Matrix(F289, 4,1, e ) ;
Jc := Evaluate(J,b) ;
s1 := Eltseq(Jc[1]) ;
s2 := Eltseq(Jc[2]) ;
s3 := Eltseq(Jc[3]) ;
s4 := Eltseq(Jc[4]) ;
d1 := [phi(s) : s in s1 ];
d2 := [phi(s) : s in s2 ];
d3 := [phi(s) : s in s3 ];
d4 := [phi(s) : s in s4 ];
s := d1 cat d2 cat d3 cat d4 ;
Jcc := Matrix(F289, 4,4, s);
Zu<[y]> := PolynomialRing(F289, 4) ;
B := Matrix(Zu, 4,1, [y[1],y[2],y[3],y[4]]) ;
JCc := RMatrixSpace(Zu, 4,4) ! Jcc ;
eqns := Eltseq( (JCc * B ) + Y ) ;
Zzz<[y]> := PolynomialRing(F289,4);
g1 := Zzz ! eqns[1] ;
g2 := Zzz ! eqns[2] ;
g3 := Zzz ! eqns[3] ;
g4 := Zzz ! eqns[4] ;
EQN := [g1,g2,g3,g4 ] ;
S1 := Scheme(AffineSpace(Zzz) , EQN) ;
PT := Points(S1) ;
PPT := [ Eltseq(a) : a in PT ] ;
PTT := [ s@@phi : s in PPT ] ;
NPT := [ [b[1] + (17^n)*s[1], b[2] + (17^n)*s[2], b[3] + (17^n)*s[3], b[4] + (17^n)*s[4] ] : s in PTT ] ;
return NPT[1] ;
end function ;

Approx := [ [lift(1,A[i])] : i in [1..28] ] ;
for j in [1..28] do ;
for i in [2..100] do ;
Approx[j][i] := lift(i, Approx[j][i-1] );
end for;
end for ;

AA := [ a[100] : a in Approx];

// we use the Hensel approximations above to compute candidates for the minimal polynomials of the coefficients a_1 and a_2 .

// we'll use the following function to compute these candidates 
// the input is an approximation x to the coefficient modulo p^k, d is the predicted degree of the polynomial 

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


MM := [ apr(a[1],100,17,6) : a in AA ];

MM := [ t : t in MM | #Factorization(t) eq 1 ];
assert #MM eq 28 ;

MM := [ t : t in MM  | Factorization(t)[1,2] eq 1 ] ;
assert #MM eq 28 ;

MM := Set(MM) ;
MM := SetToSequence(MM) ;

// we can verify that all these polynomials split over the number field  K define below. 

T<u> := PolynomialRing(Rationals());
g := -u^6 - 14*u^5 + 5*u^4 + 20*u^3 - 95*u^2 - 134*u + 139 ;
K := SplittingField(g);

RM := [] ;
for i in [1..#MM] do ;
r := MM[i] ;
t := #Roots(r,K) ;
d := Degree(r) ;
n := d -t;
RM[i] := n ;
end for ;


tt := [ 0 : i in [1..#MM] ] ;

assert RM eq tt ;

// we now verify that this is the field of definition  of the bitangents 

SK := BaseChange(S,K) ;
PTS := Points(SK) ;


assert #PTS eq 28 ;

print "The field of definition of the bitangents to X0(75)/w25 is:";
K ;


// Note that the above set PTS completely describes the bitangents to the curve and can be used to compute the 2-torsion subgroup 
// In this example, for presentation puposes, for a pair of coefficients (a1,a2)  we give a pair  (m1, r) where  m1 is the minimal polynomial (or -1* the minimal polynomial)  of a1 and a2 = r(a1). We'll search for such a relation r using lattice reduction, similar to the way we search for the minimal polynomial.
// This step can be skipped and the code in JTwoTorsion.m can be adapted to use the set PTS above, the reader should look at the code used to compute the two torsion subgroup of J0(54) for an example of this. 

mm := #MM ;
ind := [] ;

for k in [1..mm ] do ;
mi := Min([i : i in [1..28] | apr(AA[i][1],100,17,6) eq MM[k]]);
ind[k] := mi ;
end for ;




 
// we compute the rational bitangents 

MM1 := [ a : a in MM | Degree(a) eq 1 ];
ind11 := [ [ i : i in [1..#MM] | MM[i] eq f] : f in MM1] ;
ind11 := [ a[1] : a in ind11] ;
ind1 := [ ind[i] : i in  ind11] ;

MM12 := [ apr(AA[i][2],100, 17,1) : i in ind1 ];

MM1 := [ [MM1[i], MM12[i]] : i in [1..#MM1] ] ;

R1 := [ [ Roots(a[1])[1,1], Roots(a[2])[1,1]] : a in MM1 ];


// we check the above correspond to points on our scheme 

MM13 := [ apr(AA[i][3],100, 17,1) : i in ind1 ];
MM14 := [ apr(AA[i][4],100, 17,1) : i in ind1 ];

Mr := [ MM1[i] cat [ MM13[i], MM14[i] ]  : i in [1..#MM1] ] ;
RR := [ [ Roots(a, Rationals())[1,1] : a in b ] : b in Mr] ;
Ea := [ [Evaluate(e,a) : e in eqns] : a in RR];

tt := [ [0,0,0,0] : i in [1..#RR] ] ;
assert Ea eq tt ;
n1 := #RR;
print "The coefficients (a1, a2)  of the rational bitangents x -a1y - a2z  to the curve are:"; 
R1;




// we compute the bitangents defined over a quadratic field 
MM2 := [ a : a in MM | Degree(a) eq 2 ];
ind22 := [ [ i : i in [1..#MM] | MM[i] eq f] : f in MM2] ;
ind22 := [ a[1] : a in ind22] ;
ind2 := [ ind[i] : i in  ind22] ;


// we compute the relation refered to above using the following function, the input is a pair x = (x1,x2) of 17-adic approximations , a = the precision of these approximations, n the degree of the supposed relation 
// the output is a n+1 tuple of inteers (a0, a1,..,an) such that an*x2 + an-1*x1^n-1 + ... + a1x1 + a0  =0

rel := function(x,a,n) ;
ZZ := FreeAbelianGroup(n+2);
Z1 := FreeAbelianGroup(2);
Z1s := sub<Z1 |[(17^a)*Z1.1,(17^a)*Z1.2]>;
Q,pi := quo<Z1 | Z1s>;
R,r := quo<OK | P^a> ;
xx := [ s@@r : s in x ];
y := [ xx[1]^(n-t) : t in [0..n-1]] cat [xx[2]] ;
j := [ Z1 ! Eltseq(s) : s in y ];
c := [ pi(s) : s in j ] ;
d := c cat [ Q ! [1,0] ];
phi := hom< ZZ -> Q | d > ;
K := Kernel(phi) ;
G := { Eltseq(ZZ ! g) : g in Generators(K) };
W := StandardLattice(n+2);
L := sub<W | [ W ! g : g in G ] > ;
i := Index(W, L );
B := (1/1000)*((i)^(1/(n+2))) ;
b := BestApproximation(B,10^10);
v := ShortVectors(L,b) ;
V := [ Eltseq(a[1]) : a in v ]; 
Z<[z]> := PolynomialRing(Integers(),2);
Vu := [ [t[s+1]*(z[1]^(n-s)) : s in [0..n-1] ] cat [t[n+1]*z[2],t[n+2] ] : t in V ];
Rel := [ &+[a : a in b] : b in Vu ];
return  Rel;             
end function ;

Za<a> := PolynomialRing(Rationals());

R1 :=[];
R2 := [] ;
R3 := [] ;

 for i in [1..#ind2] do ;
j := ind2[i] ;
a1 := AA[j] ;
s1 := rel([a1[1],a1[2]],10,1)[1];
c := Coefficients(s1);
r := (-1/c[2])*(c[1]*a + c[3]) ;
R1[i] := r ;
s2 := rel([a1[1],a1[3]],10,1)[1];
c := Coefficients(s2);
r := (-1/c[2])*(c[1]*a + c[3]) ;
R2[i] := r ;
s3 := rel([a1[1],a1[4]],10,1)[1];
c := Coefficients(s3);
r := (-1/c[2])*(c[1]*a + c[3]) ;
R3[i] := r ;
end for ;

// we verify that the above minimal polynomials and relations do in fact define bitangents by checking that the corresponding points are solutions to the system of equations defined above  

for i in [1..#ind2] do  ;
m := MM2[i] ;
K := NumberField(m) ;
r := Roots(m,K) ;
r := [a[1] : a in r] ;
r1 := R1[i] ;
r2 := R2[i];
r3 := R3[i];
pts := [ [a, Evaluate(r1, a), Evaluate(r2, a), Evaluate(r3,a)] : a in r];
t := [ [ Evaluate(e,a) : e in eqns] : a in pts ];
for j in [1..#t] do ;
assert t[j] eq [0,0,0,0];
end for ;
end for ;

print "The bitangents defined over quadratic fields are defined by equations  of the form x = a_1y + a_2z, and below we give (m1,r), m1  the minimal polynomial of a1 and r a relation such that a_2 = r(a_1).";

for i in [1..#ind2] do ;
<MM2[i], R1[i]> ;
end for ;

n2 := 2*(#ind2);
// we repeat the above to compute the bitangents defined over cubic fields 

MM3 := [ a : a in MM | Degree(a) eq 3 ];
assert #MM3 eq 0 ;
print "There are no bitangents defined over cubic fields."; 

// we repeat the above to find the bitangent defined over a number field of degree 6.

MM3 := [ a : a in MM | Degree(a) eq 6 ];
ind33 := [ [ i : i in [1..#MM] | MM[i] eq f] : f in MM3] ;
ind33 := [ a[1] : a in ind33] ;
ind3 := [ ind[i] : i in  ind33] ;

n3 := 6*(#ind3);


Za<a> := PolynomialRing(Rationals());


R1 :=[];
R2 := [] ;
R3 := [] ;

 for i in [1..#ind3] do ;
j := ind3[i] ;
a1 := AA[j] ;
s1 := rel([a1[1],a1[2]],52,5)[1];
c := Coefficients(s1);
r := (-1/c[6])*(c[1]*a^5 + c[2]*a^4 + c[3]*a^3 + c[4]*a^2 + c[5]*a + c[7]) ;
R1[i] := r ;
s2 := rel([a1[1],a1[3]],52,5)[1];
c := Coefficients(s2);
r := (-1/c[6])*(c[1]*a^5 + c[2]*a^4 + c[3]*a^3 + c[4]*a^2 + c[5]*a + c[7]) ;
R2[i] := r ;
s3 := rel([a1[1],a1[4]],52,5)[1];
c := Coefficients(s3);
r := (-1/c[6])*(c[1]*a^5 + c[2]*a^4 + c[3]*a^3 + c[4]*a^2 + c[5]*a + c[7]) ;
R3[i] := r ;
end for ;

for i in [1..#ind3] do  ;
m := MM3[i] ;
K := NumberField(m) ;
r := Roots(m,K) ;
r := [a[1] : a in r] ;
r1 := R1[i] ;
r2 := R2[i];
r3 := R3[i];
pts := [ [a, Evaluate(r1, a), Evaluate(r2, a), Evaluate(r3,a)] : a in r];
t := [ [ Evaluate(e,a) : e in eqns] : a in pts ];
for j in [1..#t] do ;
assert t[j] eq [0,0,0,0];
end for ;
end for ;

print "The bitangents defined over a number  field of degree 6  are defined by equations  of the form x = a_1y + a_2z, and below we state (m1,r), the minimal polynomil of a1 and r, a relation such that a_2 = r(a_1)." ;

MM3;
R1;
for i in [1..#ind3] do ;
<MM3[i], R1[i]> ;
end for ;

assert n1 + n2 + n3 eq 28 ;
print "This is the complete list of bitangents to the curve.";

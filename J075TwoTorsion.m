// This is the MAGMA code used to compute the 2-torsion subgroup of J0(75) 

Z<U> := PolynomialRing(Integers());

f := U^12 + 6*U^11 - 91*U^10 - 90*U^9 + 3370*U^8 - 14834*U^7 + 38281*U^6 - 79946*U^5 + 123730*U^4 - 130890*U^3 + 127789*U^2 - 57586*U + 58381;
K := NumberField(f) ;

// all the quadritangents found are defined over this number field 

PP := PrimesUpTo(150,K) ;
P := PP[19] ;

// this is a prime ideal of norm 121
OK := Integers(K) ;


RR<[x]> := PolynomialRing(Rationals(),5) ;
f1 := x[1]*x[3] - x[2]^2 - x[4]^2 - 4*x[5]^2 ;
f2 := x[1]*x[4] - x[2]*x[3]  + x[2]*x[5] + x[4]*x[5] - 3*x[5]^2 ;
f3 := x[1]*x[5] - x[2]*x[4] - x[3]*x[5] ;
P4 := ProjectiveSpace(RR) ;

X := Curve(P4, [f1,f2,f3] ) ;

// we reduce our curve modulo P, and we'll use the fact that the induced map on the Jacobian is injective on torsion 

F121,pi := ResidueClassField(P) ;
X121 := ChangeRing(X, F121) ;
Cl121, phi, psi := ClassGroup(X121) ;
Z := FreeAbelianGroup(1) ;
degr := hom< Cl121 -> Z | [ Degree(phi(g)) : g in OrderedGenerators(Cl121) ] >;


J121 := Kernel(degr) ; 
// this is the set of points of the Jacobian over the finite field F_121

// we now define the equations of the quadritangent planes computed using X075Planes.m 

f1 := f;
R := Roots(f1,K) ;
R1 := R ;
Z<u> := PolynomialRing(K);
S1 := (-1/141639721686689040)*(1622551063058*u^11 + 12639757545665*u^10 - 124528372282287*u^9 - 365357454364857*u^8 + 4773709008942128*u^7 - 15606572136946403*u^6 + 35791517616953403*u^5 - 70223542383995731*u^4 + 91562316359380564*u^3 - 78677279262041201*u^2 + 77253656939824629*u - 
10294852788765071 ) ;
S2 := (-1/70819860843344520)*(1302412942954*u^11 + 9339108581185*u^10 - 106104169871511*u^9 - 230508625151001*u^8 + 3990547234584040*u^7 - 15061292130592219*u^6 + 36712870645335075*u^5 - 75233014631294963*u^4 + 105419530746350132*u^3 -  83119161971186305*u^2 + 103496216790983997*u + 
19214821033550825);
S3 := (1/694312361209260)*(28676117706*u^11 + 215479079675*u^10 - 2261103354449*u^9 - 5841824308979*u^8 + 85924080818884*u^7 - 300665335956261*u^6 + 710827335904789*u^5 - 1426044676620497*u^4 + 1931194579467948*u^3 - 1586239619933603*u^2 + 1077745224386903*u - 1301174053351993);
coeff1 := [ [-1, R[i,1], Evaluate(S1, R[i,1]), Evaluate(S2,R[i,1]),  Evaluate(S3, R[i,1]) ] : i in [1..12] ];

Z<U> := PolynomialRing(Integers());
f1 := U^12 + U^11 + 9*U^10 + 40*U^9 + 125*U^8 + 431*U^7 + 1456*U^6 + 639*U^5 + 3895*U^4 + 1200*U^3 + 2859*U^2 + 9*U + 711 ;

R := Roots(f1, K) ;   
R2 := R ;                                                       

Z<u> := PolynomialRing(K); 

S1 := (1/202093678756440)*(31170868775*u^11 + 2980773224*u^10 + 216795559605*u^9 + 957429111941*u^8 + 
    2624059395016*u^7 + 8123871102031*u^6 + 28621146996425*u^5 - 
    27997400094900*u^4 + 61774433256743*u^3 - 101497172571957*u^2 - 
    30299240860572*u- 65505228496371);
S2 := (1/67364559585480)*(40087728559*u^11 + 58541592208*u^10 + 327286175019*u^9 + 1687350552805*u^8 + 5438483245694*u^7 + 17603054304743*u^6 + 59289585628945*u^5 + 32117274981912*u^4 + 99073013306845*u^3 + 109752651751827*u^2 + 39338993277258*u + 23930256987153);
S3 := (-1/202093678756440)*(459044962673*u^11 + 126506149952*u^10 + 3740752976877*u^9 + 15714131893187*u^8 +
    44107882097254*u^7 + 157655613428317*u^6 + 533596098579095*u^5 - 
    148827339389304*u^4 + 1646864492918651*u^3 - 259720034066643*u^2 + 
    723559950639474*u - 262869802509213);

coeff2 := [ [-1, R[i,1], Evaluate(S1, R[i,1]), Evaluate(S2,R[i,1]),  Evaluate(S3, R[i,1]) ] : i in [1..12] ]; 

Z<U> := PolynomialRing(Integers());
f1 :=  U^12 - 2*U^11 + 21*U^10 - 50*U^9 + 130*U^8 - 1002*U^7 + 2169*U^6 - 7362*U^5 + 16650*U^4 - 21330*U^3 + 50301*U^2 - 22842*U + 53541;
R := Roots(f1, K) ;       
R3 := R ;
Z<u> := PolynomialRing(K); 
S1 := (-1/989352072720)*(31780654*u^11 + 35184793*u^10 + 487829205*u^9 + 92792059*u^8 - 174066092*u^7 - 
    24459595611*u^6 - 20346101625*u^5 - 64385877375*u^4 + 113713491096*u^3 + 
    646293110247*u^2 + 738498973977*u + 1579202221077);
S2 := (1/14990182920)*(1112348*u^11 - 2828753*u^10 + 18039343*u^9 - 67448645*u^8 + 66156922*u^7 - 
    1115875817*u^6 + 2484057557*u^5 - 4549266183*u^4 + 18404993250*u^3 - 
    11014499547*u^2 + 28019583279*u- 27505059759);
S3 := (-1/247338018180)*(2463415*u^11 - 64266821*u^10 + 53734557*u^9 - 1159298672*u^8 + 1178622259*u^7 - 
    6182153175*u^6 + 51160000503*u^5 - 42869953332*u^4 + 246825643077*u^3 - 
    504885797649*u^2 + 340411655295*u - 748758560202);


coeff3 := [ [-1, R[i,1], Evaluate(S1, R[i,1]), Evaluate(S2,R[i,1]),  Evaluate(S3, R[i,1]) ] : i in [1..12] ];

Z<U> := PolynomialRing(Integers());
f1 := U^12 - 2*U^11 + 21*U^10 + 30*U^9 + 250*U^8 + 238*U^7 + 1369*U^6 - 602*U^5 + 5770*U^4 - 8250*U^3 + 12261*U^2 - 20042*U + 10981;
R := Roots(f1, K) ;      
R4 := R ;
S1 := (-1/157202076240)*(9745328*u^11 - 17139225*u^10 + 239277331*u^9 + 279059519*u^8 + 3183501206*u^7 + 
    4859289751*u^6 + 22254727097*u^5 + 6337541909*u^4 + 96208059166*u^3 - 
    57723043371*u^2 + 236633195771*u - 316718987547 ) ;
S2 := (1/39300519060)*(9668974*u^11 - 686748*u^10 + 164436374*u^9 + 675526051*u^8 + 3075298126*u^7 + 
    6622159106*u^6 + 18558678826*u^5 + 21591154141*u^4 + 58285309094*u^3 + 
    31275315354*u^2 + 31646937670*u - 55390969533 );

S3 := (-1/78601038120)*(21427516*u^11 - 49380907*u^10 + 464907767*u^9 + 612594113*u^8 + 4957582866*u^7 +
    5885083309*u^6 + 32453017757*u^5 + 3744433923*u^4 + 133521145042*u^3 - 
    127267891241*u^2 + 208754169895*u - 157211951845);
coeff4 := [ [-1, R[i,1], Evaluate(S1, R[i,1]), Evaluate(S2,R[i,1]),  Evaluate(S3, R[i,1]) ] : i in [1..12] ]; 
coeff := coeff1  cat coeff2 cat coeff3 cat coeff4 ;
n := #coeff;

// we now clear denominators and scale appropriately to ensure that all our equations are defined over the ring of integers OK 

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

FX121 := FunctionField(X121) ;
CQT := [ [MonomialCoefficient(a,X[1]),  MonomialCoefficient(a,X[2]) , MonomialCoefficient(a,X[3]) , MonomialCoefficient(a,X[4]), 
MonomialCoefficient(a,X[5])] : a in QT ] ;
F121QT := [ [ pi(s) : s in a ] : a in CQT ] ;

ZU<[U]> := PolynomialRing(F121, 5) ;

F121QT := [ a[1]*U[1] + a[2]*U[2] + a[3]*U[3] + a[4]*U[4] + a[5]*U[5] : a in F121QT];
F121QTT := [ Evaluate(a,[ FX121.1, FX121.2, FX121.3, FX121.4, 1]) : a in F121QT] ;
F121QTT := [ a : a in F121QTT | a ne 0 ] ;

nn := #F121QTT ;
D := [Divisor(F121QTT[i]/F121QTT[1]) : i in [2..nn] ]  ;

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
hh := hom< ZN -> J121 | [ a : a in H ] > ;
ihh := Image(hh) ;

print "The 2-torsion subgroup generated by the computed quadritangets is:";
ihh ;

assert #ihh eq 2^10;
print "Thus we have computed the entire 2-torsion subgroup.";

// we take Galois invariants to obtain the rational part 


// we find that the Galois group of K has 2-generatos, which act as  cpt1 and cpt2 below on our generators 

cpt1 := [ ZN.4 - ZN.10, ZN.6 - ZN.10, ZN.11 - ZN.10 , ZN.9 - ZN.10 , ZN.3 - ZN.10, ZN.8 - ZN.10, - ZN.10, ZN.5 - ZN.10, ZN.7 - ZN.10, ZN.1 - ZN.10, ZN.2 - ZN.10 , ZN.15 - ZN.10, ZN.23 - ZN.10, ZN.13 - ZN.10, ZN.16 - ZN.10, ZN.17 - ZN.10, ZN.22 - ZN.10, ZN.20 - ZN.10, ZN.12 - ZN.10, ZN.14 - ZN.10, ZN.18 - ZN.10, ZN.19 - ZN.10, ZN.21 - ZN.10, ZN.31 - ZN.10, ZN.30 - ZN.10, ZN.29 - ZN.10, ZN.25 - ZN.10, ZN.32 - ZN.10, ZN.33 - ZN.10, ZN.28 - ZN.10, ZN.35 - ZN.10, ZN.34 - ZN.10, ZN.24 - ZN.10, ZN.27 - ZN.10, ZN.26 - ZN.10, ZN.43 - ZN.10, ZN.42 - ZN.10, ZN.36 - ZN.10, ZN.38 - ZN.10, ZN.41 - ZN.10, ZN.37 - ZN.10, ZN.44 - ZN.10, ZN.47 - ZN.10, ZN.45 - ZN.10, ZN.40 - ZN.10, ZN.39 - ZN.10, ZN.46 - ZN.10];

conj1 := hom< ZN -> ZN | cpt1>;
mu := hom< ZN -> J121 | [ hh(ZN.i) - hh(conj1(ZN.i)) : i in [1..47]]>;
ker1 := Kernel(mu);
imKer1 := sub<J121 | [hh(k) : k in Generators(ker1)]>;

cpt2 := [ ZN.8 - ZN.3, ZN.9 - ZN.3, - ZN.3, ZN.6 - ZN.3, ZN.10 - ZN.3, ZN.4 - ZN.3, ZN.11 - ZN.3, ZN.1 - ZN.3, ZN.2 - ZN.3, ZN.5 - ZN.3, ZN.7 - ZN.3, ZN.23 - ZN.3, ZN.15 - ZN.3, ZN.16 - ZN.3, ZN.13 - ZN.3, ZN.14 - ZN.3, ZN.20 - ZN.3, ZN.22 - ZN.3, ZN.21 - ZN.3, ZN.17 - ZN.3, ZN.19 - ZN.3, ZN.18 - ZN.3, ZN.12 - ZN.3 ,ZN.34 - ZN.3, ZN.29 - ZN.3, ZN.30 - ZN.3, ZN.33 - ZN.3, ZN.35 - ZN.3, ZN.25 - ZN.3, ZN.26 - ZN.3, ZN.32 - ZN.3, ZN.31 - ZN.3, ZN.27 - ZN.3, ZN.24 - ZN.3, ZN.28 - ZN.3, ZN.37 - ZN.3, ZN.36 - ZN.3, ZN.42 - ZN.3, ZN.44 - ZN.3, ZN.47 - ZN.3, ZN.43 - ZN.3, ZN.38 - ZN.3, ZN.41 - ZN.3, ZN.39 - ZN.3, ZN.46 - ZN.3, ZN.45 - ZN.3, ZN.40 - ZN.3 ];

conj2 := hom< ZN -> ZN | cpt2>;
mu := hom< ZN -> J121 | [ hh(ZN.i) - hh(conj2(ZN.i)) : i in [1..47]]>;
ker2 := Kernel(mu);
imKer2 := sub<J121 | [hh(k) : k in Generators(ker2)]>;

rat2tors := imKer1 meet imKer2 ;

print "The rational 2-torsion subgroup is:";
rat2tors ;


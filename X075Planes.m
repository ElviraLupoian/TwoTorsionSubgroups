// This is the MAGMA code used to compute the quadritangents to X0(75) used in our computation of the 2-torsion subgroup of its Jacobian 

// the coefficinets of the equations defining the quadritangents used are solutions of the system of equations:

Zx<[x]> := FunctionField(Integers(), 5) ;
f1 := x[1]*x[3] - x[2]^2 - x[4]^2 - 4*x[5]^2 ;
f2 := x[1]*x[4] - x[2]*x[3]  + x[2]*x[5] + x[4]*x[5] - 3*x[5]^2 ;
f3 := x[1]*x[5] - x[2]*x[4] - x[3]*x[5] ;

Za<[a]> := PolynomialRing(Integers(),19) ;
ZX<[X]> := FunctionField(Za,4) ;
F1 := Evaluate(f1, X cat [1]) ;
F2 := Evaluate(f2, X cat [1]) ;
F3 := Evaluate(f3, X cat [1]) ;

G1 := Evaluate(F1, [ a[1]*X[2] + a[2]*X[3] + a[3]*X[4] + a[4], X[2], X[3], X[4] ] ) ;
G2 := Evaluate(F2, [ a[1]*X[2] + a[2]*X[3] + a[3]*X[4] + a[4], X[2], X[3], X[4] ] ) ;
G3 := Evaluate(F3, [ a[1]*X[2] + a[2]*X[3] + a[3]*X[4] + a[4], X[2], X[3], X[4] ] ) ;

X33 := (-X[2]*X[4] + a[1]*X[2] + a[3]*X[4] + a[4]) ;
h := 1 - a[2] ;
h := ZX ! h ;
h := 1/h ;

X3 := h*X33 ;

g1 := Evaluate(G1, [ 1, X[2], X3, X[4]] ) ;
g2 := Evaluate(G2, [ 1, X[2], X3, X[4]] ) ;
gg1 := Numerator(g1) ;
gg2 := Numerator(g2) ;

Zxy<x,y> := FunctionField(Za, 2) ;

e1 := Evaluate(gg1, [ 1, x, 1, y] ) ;
e2 := Evaluate(gg2, [ 1, x, 1, y] ) ;

t1 := a[2]*x^2 + (-a[2]*a[3] - a[3])*x + (-a[2]^2 + 2*a[2] + a[3]^2 - 1) ;
t2 := (-a[1]*a[2] - a[1])*x^2 + (2*a[1]*a[3] - a[2]*a[4] - a[4])*x + 2*a[3]*a[4] ;

t3 := (a[1]^2 - a[2]^2 + 2*a[2] - 1)*x^2 + 2*a[1]*a[4]*x - 4*a[2]^2 + 8*a[2] + a[4]^2 - 4 ;


s1 := a[2]*x - a[3] ;
s2 := -x^2 + (-a[1] + a[3])*x + (a[2] - a[4] - 1) ;
s3 :=  a[1]*x^2 + (a[2] + a[4] - 1)*x - 3*a[2] + 3 ;

alpha := t1*s2 - s1*t2 ;

beta := t1*s3 -s1*t3 ;

Y := -beta/alpha ;

H := Evaluate(e2, [x, Y] ) ;
nH := Numerator(H) ;
Zt<t> := PolynomialRing(Za);
tHz := Evaluate(nH, [ t,1]) ;


H := Factorization(tHz)[3,1] ;
h := t^4 + a[5]*t^3 + a[6]*t^2 + a[7]*t + a[8] ;
l := MonomialCoefficient(H,t^8) ;
E := H - l*h^2 ;
Eq := Coefficients(E) ;
Eq[9] := a[9]*(Discriminant(h)) + 1 ;

Eq[10] := a[10]*(a[2] -1 ) + 1 ;
Eq[11] := a[11]*a[2] + 1 ;
Eq[12] := a[12]*a[3] + 1 ;
Eq[13] := a[13]*(a[2] + 1) + 1 ;
Eq[14] := a[14]*(a[2] - a[3] -1 ) + 1 ;
Eq[15] := a[15]*(a[2] + a[3] - 1) + 1 ;
Eq[16] := a[16]*(a[1]*a[2]^2 + 2*a[2]*a[3] + a[3]) + 1 ;
Eq[17] := a[17]*(2*a[1]*a[2]*a[3] - a[2]^2*a[4] - 2*a[2]^2 + a[2]*a[3]^2 + 3*a[2] + 2*a[3]^2- 1) + 1;
Eq[18] := a[18]*(a[1]*a[2]^2 - 2*a[1]*a[2] + a[1]*a[3]^2 + a[1] - 2*a[2]^2*a[3] - 2*a[2]*a[3]*a[4] + 2*a[2]*a[3] + a[3]^3) + 1 ;
Eq[19] := a[19]*( a[2]^3 - a[2]^2*a[4] - 3*a[2]^2 - a[2]*a[3]^2 + 2*a[2]*a[4] + 3*a[2] -  a[3]^2*a[4] + a[3]^2 - a[4] - 1) + 1 ;

CC := ComplexField(2000) ;
ZZa<[a]> := PolynomialRing(CC,19) ;
EE := [] ;
for i in [1..19] do ;
EE[i] := ZZa ! Eq[i] ;
end for ;

// EE is the list of equations defining our zero-dimensional scheme 


J := JacobianMatrix(EE) ;
im := CC.1 ;
e := Exp(1) ;
e := CC ! e ;

load "general.m";

// the first approximation we use:

a1 := [ 3.66300314429997260000000000000 + 0.562539534077005300000000000000*im, 
-2.66156941367538870000000000000 - 2.54969885172509430000000000000*im, 
-1.66921191810016260000000000000 + 2.74291918955430520000000000000*im, 
1.32934760115097420000000000000 + 1.79393897981887850000000000000*im, 
-7.76534990627579300000000000000 + 1.33327118995548590000000000000*im, 
-4.22581752801258050000000000000 - 11.1585974322272940000000000000*im, 
-31.4423519940689040000000000000 + 2.63021186721544130000000000000*im, 
-24.1057159135274500000000000000 - 56.5473569180709500000000000000*im, 
-0.000467589527260020993713435392557 - 0.000294760220770526234099015423787*im, 
0.183924017158922100000000000000 - 0.128073730789675100000000000000*im, 
0.195920933035342410000000000000 - 0.187686022924084700000000000000*im, 
0.161904199859488060000000000000 + 0.266047786891836300000000000000*im, 
0.179400711828289570000000000000 - 0.275292615031620600000000000000*im, 
0.0622975883350449300000000000000 - 0.165491053027640030000000000000*im, 
0.187343632705989200000000000000 + 0.00679048675019389500000000000000*im, 
-0.00646975704299248160000000000000 + 0.0192386380701610320000000000000*im, 
-0.00945154955157480300000000000000 - 0.00531039094352075100000000000000*im, 
-0.00893457889123702700000000000000 + 0.00395809912210298450000000000000*im, 
-0.00178998267489239880000000000000 - 0.00623743951336383500000000000000*im  ] ;

a1 := [ CC ! b : b in a1];
N := [a1] ;
J := JacobianMatrix(EE) ;


// we perform 700 steps of Newton-Raphson 


for i in [2..700] do ;
a := N[i-1] ;
a := [ CC ! b : b in a] ;
N[i] := nr(a,19) ;
end for ;

pt := N[700] ;

// we search for the minimal polynomial of a[1] and for relations for a[2], a[3] and a[4] in terms of a[1]

m1 := minpoly(pt[1],700,12) ;
r1 := re([pt[1],pt[2]], 700, 12) ;
r2 := re([pt[1], pt[3]], 700, 12) ;
r3 := re([pt[1], pt[4]], 700, 12) ;

// we verify that roots of the above equation are points on our scheme 

SK := SplittingField(m1) ;
R := Roots(m1,SK);
R := [ a[1] : a in R ];
r := [ r1,r2, r3] cat [ re([pt[1],pt[i]],700,12) : i in [5..19] ] ;
PT := [ [a] cat [ Evaluate(s,a) : s in r ] : a in R ] ;

ZX<[X]> := PolynomialRing(Integers(),19) ;
TT := [ Evaluate(e, X) : e in EE];  // we base change the field of definition of the equations in EE 
EPT := [ [ Evaluate(e,a) : e in TT ] : a in PT] ;
NET := [ a : a in EPT | a ne [ 0 : i in [1..19] ]] ;

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


// the second orbit of quadritangents corresponds to the following approximation 


a1 := [ 0.430900702718914300000000000000 + 1.46104950105113100000000000000*im, 
0.401649919855882500000000000000 - 0.243530565905747930000000000000*im, 
0.482727146937039400000000000000 + 1.38776470243943440000000000000*im, 
2.58797964536951630000000000000 - 0.653025495563826400000000000000*im, 
0.148836974546938280000000000000 - 3.18159026154313300000000000000*im, 
2.14041752031117750000000000000 + 1.95128858549174120000000000000*im, 
9.29476497750391100000000000000 - 14.7269467191707460000000000000*im, 
-11.1163420650247340000000000000 - 9.35899623201824600000000000000*im, 
-0.00159202191837041150000000000000 - 0.00237811366022264000000000000000*im, 
1.43375780497663530000000000000 - 0.583544418567879200000000000000*im, 
-1.82046998906895260000000000000 - 1.10379727403274300000000000000*im, 
-0.223596964733433930000000000000 + 0.642806142555985300000000000000*im, 
-0.692538859088083900000000000000 - 0.120325609038513070000000000000*im, 
0.282276485882383900000000000000 - 0.425942092063391840000000000000*im, 
0.0874182667437668800000000000000 + 0.865113539246400900000000000000*im, 
-0.209415190821907650000000000000 + 0.260310231679492400000000000000*im, 
0.0843976274206999400000000000000 + 0.104056611316232730000000000000*im, 
0.0972670209398633900000000000000 - 0.0435264756397054850000000000000*im, 
-0.0419868321615776400000000000000 - 0.191188937333945620000000000000*im ];


a1 := [ CC ! b : b in a1];
N := [a1] ;
J := JacobianMatrix(EE) ;


// we perform 700 steps of Newton-Raphson 


for i in [2..700] do ;
a := N[i-1] ;
a := [ CC ! b : b in a] ;
N[i] := nr(a,19) ;
end for ;

pt := N[700] ;

// we search for the minimal polynomial of a[1] and for relations for a[2], a[3] and a[4] in terms of a[1] 

m1 := minpoly(pt[1],700,12) ;
r1 := re([pt[1],pt[2]], 700, 12) ;
r2 := re([pt[1], pt[3]], 700, 12) ;
r3 := re([pt[1], pt[4]], 700, 12) ;

// we verify that roots of the above equation are points on our scheme 

SK := SplittingField(m1) ;
R := Roots(m1,SK);
R := [ a[1] : a in R ];
r := [ r1,r2, r3] cat [ re([pt[1],pt[i]],700,12) : i in [5..19] ] ;
PT := [ [a] cat [ Evaluate(s,a) : s in r ] : a in R ] ;
ZX<[X]> := PolynomialRing(Integers(),19) ;
TT := [ Evaluate(e, X) : e in EE];  // we base change the field of definition of EE
EPT := [ [ Evaluate(e,a) : e in TT ] : a in PT] ;
NET := [ a : a in EPT | a ne [ 0 : i in [1..19] ]] ;

assert NET eq [] ;

print "One orbit of 4-tangent planes of the form -x[1] + a[1]x[2] + 
a[2]x[3] + a[3]x[4] + a[4]x[5] is as follows:";

print " The minimal polynomial of a[1] is:" ;
m1;

print " and relations for a[2], a[3], a[4] resp. in terms of a[1] are:";
r1;
r2;
r3;

O2 := [ m1, r1, r2, r3];

// the third orbit of quadritangents corresponds to the following approximation 

a1 := [ 0.6185589238559607 + 4.256034574480145*CC.1, 11.89647267150720 - 
9.128437813098847*CC.1, -19.31799642601400 - 16.78412800596857*CC.1, 
-7.093507840856364 + 30.78496905768612*CC.1, -23.47871830013747 - 
5.361752466409918*CC.1, 152.7638798372919 + 5.262679763159169*CC.1, 
-182.1535761617178 + 377.4126773086385*CC.1, -245.5434021631166 - 
149.1026592779166*CC.1, 2.988908738616839E-5 + 1.886253298429669E-5*CC.1, 
-0.05392651750707679 - 0.04517653339576050*CC.1, -0.05290743971094004 - 
0.04059709853394759*CC.1, 0.02949799085990116 - 0.02562885112892778*CC.1, 
-0.05165875344131205 - 0.03656532528720598*CC.1, -0.03110008330584073 + 
0.007880085597184602*CC.1, 0.01134390507794949 - 0.03490457256269228*CC.1, 
-0.005284523847554025 + 0.001506933109471823*CC.1, -0.0001345961911046843 + 
0.0004712645750597950*CC.1, -0.06394962904934870 + 0.02452953552801368*CC.1, 
-0.03314217505608918 - 0.04893789729657522*CC.1 ];


a1 := [ CC ! b : b in a1];
N := [a1] ;
J := JacobianMatrix(EE) ;


// we perform 700 steps of Newton-Raphson 


for i in [2..700] do ;
a := N[i-1] ;
a := [ CC ! b : b in a] ;
N[i] := nr(a,19) ;
end for ;

pt := N[700] ;


// we search for the minimal polynomial of a[1] and for relations for a[2], a[3] and a[4] in terms of a[1] 

m1 := minpoly(pt[1],700,12) ;
r1 := re([pt[1],pt[2]], 700, 12) ;
r2 := re([pt[1], pt[3]], 700, 12) ;
r3 := re([pt[1], pt[4]], 700, 12) ;

// we verify that roots of the above equation are points on our scheme 

SK := SplittingField(m1) ;
R := Roots(m1,SK);
R := [ a[1] : a in R ];
r := [ r1,r2, r3] cat [ re([pt[1],pt[i]],700,12) : i in [5..19] ] ;
PT := [ [a] cat [ Evaluate(s,a) : s in r ] : a in R ] ;
ZX<[X]> := PolynomialRing(Integers(),19) ;
TT := [ Evaluate(e, X) : e in EE];  // we base change the field of definition of EE
EPT := [ [ Evaluate(e,a) : e in TT ] : a in PT] ;
NET := [ a : a in EPT | a ne [ 0 : i in [1..19] ]] ;

assert NET eq [] ;

print "One orbit of 4-tangent planes of the form -x[1] + a[1]x[2] + 
a[2]x[3] + a[3]x[4] + a[4]x[5] is as follows:";

print " The minimal polynomial of a[1] is:" ;
m1;

print " and relations for a[2], a[3], a[4] resp. in terms of a[1] are:";
r1;
r2;
r3;

O3 := [ m1, r1, r2, r3];





// our final orbit corresponds to the following approximation

a1 := [ -0.5882713833996805 - 2.382955706381257*CC.1, 0.8880075936507242 - 
0.8506830331837264*CC.1, 0.5983157508621996 - 0.9151491809980560*CC.1, 
0.6975004074639968 + 1.076644468357341*CC.1, -2.415194726538166 - 
1.632191966982423*CC.1, 1.365560648484030 - 0.5103535342828849*CC.1, 
-1.891854265853361 + 0.8798995160844088*CC.1, 1.982893139929113 + 
0.9197622259367405*CC.1, -0.03307651687419794 + 0.0001069090573177763*CC.1, 
0.1521214476024408 - 1.155499186749822*CC.1, -0.5872215131875489 - 
0.5625395340770060*CC.1, -0.5004817064038047 - 0.7655078828527093*CC.1, 
-0.4402762216455131 - 0.1983760621131080*CC.1, 1.396338011118258 + 
0.1267288453745246*CC.1, -0.1449688920348369 - 0.5263797070357437*CC.1, 
0.1726734049501653 - 0.1379778998216427*CC.1, 0.1148852364467978 + 
0.01943976309041539*CC.1, 0.09398467078821068 + 0.1481939069507564*CC.1, 
-0.2026321534588531 + 0.4080620911052110*CC.1 ];


a1 := [ CC ! b : b in a1];
N := [a1] ;
J := JacobianMatrix(EE) ;


// we perform 700 steps of Newton-Raphson 


for i in [2..700] do ;
a := N[i-1] ;
a := [ CC ! b : b in a] ;
N[i] := nr(a,19) ;
end for ;

pt := N[700] ;

// we search for the minimal polynomial of a[1] and for relations for a[2], a[3] and a[4] in terms of a[1] 

m1 := minpoly(pt[1],700,12) ;
r1 := re([pt[1],pt[2]], 700, 12) ;
r2 := re([pt[1], pt[3]], 700, 12) ;
r3 := re([pt[1], pt[4]], 700, 12) ;

// we verify that roots of the above equation are points on our scheme 

SK := SplittingField(m1) ;
R := Roots(m1,SK);
R := [ a[1] : a in R ];
r := [ r1,r2, r3] cat [ re([pt[1],pt[i]],700,12) : i in [5..19] ] ;
PT := [ [a] cat [ Evaluate(s,a) : s in r ] : a in R ] ;
ZX<[X]> := PolynomialRing(Integers(),19) ;
TT := [ Evaluate(e, X) : e in EE];  // we base change the field of definition of EE
EPT := [ [ Evaluate(e,a) : e in TT ] : a in PT] ;
NET := [ a : a in EPT | a ne [ 0 : i in [1..19] ]] ;
assert NET eq [] ;

print "One orbit of 4-tangent planes of the form -x[1] + a[1]x[2] + 
a[2]x[3] + a[3]x[4] + a[4]x[5] is as follows:";

print " The minimal polynomial of a[1] is:" ;
m1;

print " and relations for a[2], a[3], a[4] resp. in terms of a[1] are:";
r1;
r2;
r3;

O4 := [ m1, r1, r2, r3];


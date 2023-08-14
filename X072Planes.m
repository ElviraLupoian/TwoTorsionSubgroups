//  This is the MAGMA code used to compute the quadritangenets used in our computations of the two-torsion subgroup of the modular jacobian J0(72) 

// the coefficients of the defining equations of our quadritangents are solutions to the following system of equations 

Zx<[x]> := FunctionField(Integers(), 5) ;
f1 := x[1]*x[3] -x[2]^2 + x[4]^2 ;
f2 := x[1]*x[4] - x[3]^2 ;
f3 := x[1]*x[5] - x[3]*x[4] -2*x[5]^2 ;

Za<[a]> := PolynomialRing(Integers(),16) ;
ZX<[X]> := FunctionField(Za,4) ;
F1 := Evaluate(f1, X cat [1]) ;
F2 := Evaluate(f2, X cat [1]) ;
F3 := Evaluate(f3, X cat [1]) ;


G1 := Evaluate(F1, [ a[1]*X[2] + a[2]*X[3] + a[3]*X[4] + a[4], X[2], X[3], X[4] ] ) ;
G2 := Evaluate(F2, [ a[1]*X[2] + a[2]*X[3] + a[3]*X[4] + a[4], X[2], X[3], X[4] ] ) ;
G3 := Evaluate(F3, [ a[1]*X[2] + a[2]*X[3] + a[3]*X[4] + a[4], X[2], X[3], X[4] ] ) ;

G3 = a[1]*X[2] - X[3]*X[4] + a[2]*X[3] + a[3]*X[4] + a[4] - 2 ;

hx :=  - X[3]*X[4] + a[2]*X[3] + a[3]*X[4] + a[4] - 2;
a1 := ZX ! ( -a[1] ) ;
X2 := hx/a1 ;

g1 := Evaluate(G1, [ 1, X2, X[3], X[4]] ) ;
g2 := Evaluate(G2, [ 1, X2, X[3], X[4]] ) ;



g1 := (-X[3]^2*X[4]^2 + (a[1]^2 + 2*a[2])*X[3]^2*X[4] - a[2]^2*X[3]^2 + 
    2*a[3]*X[3]*X[4]^2 + (-2*a[2]*a[3] + 2*a[4] - 4)*X[3]*X[4] + (2*a[1]^2 - 
    2*a[2]*a[4] + 4*a[2])*X[3] + (a[1]^2 - a[3]^2)*X[4]^2 + (-2*a[3]*a[4] + 
    4*a[3])*X[4] - a[4]^2 + 4*a[4] - 4)/a[1]^2;

g2 := -X[3]^3 + X[3]*X[4]^2 + 2*X[4];

g12 := -X[3]^2 +  2*a[3]*X[3]+ (a[1]^2 - a[3]^2) ;

g11 :=  (a[1]^2 + 2*a[2])*X[3]^2 + (-2*a[2]*a[3] + 2*a[4] - 4)*X[3] + (-2*a[3]*a[4] + 4*a[3]);
g10 := - a[2]^2*X[3]^2 + (2*a[1]^2 - 2*a[2]*a[4] + 4*a[2])*X[3] - a[4]^2 + 4*a[4] - 4 ;

g22 := X[3] ;
g21 := 2 ;
g20 := -X[3]^2 ;


ng1 := Numerator(g1) ;
ng2 := Numerator(g2) ; 

 alpha := g22*g11 -g12*g21 ;
beta := g22*g10 - g12*g20 ;


X4 := -beta/alpha ;

hh1 := Evaluate(ng1, [ 1, 1, X[3], X4] ) ;
hh2 := Evaluate(ng2, [ 1, 1, X[3], X4] ) ;

h1 := Numerator(hh1) ;
H := Factorization(h1)[3,1] ;

ZT<T> := PolynomialRing(Za) ; 
H := Evaluate(H, [ 1, 1, T, 1 ] ) ;
h := T^4 + a[5]*T^3 + a[6]*T^2 + a[7]*T + a[8];


E := H - h^2 ;
CE := Coefficients(E) ;
d := Discriminant(h) ;
e9 := a[9]*d + 1 ;
CE[2] := Factorization(CE[2])[2,1];
CE[4] := Factorization(CE[4])[2,1];
CE[6] := Factorization(CE[6])[2,1];

CE[9] := a[9]*d + 1;
CE[10] := a[10]*a[1] + 1 ;



CE[11] := a[11]*a[3] + 1 ;
CE[12] := a[12]*(a[1] - a[3]) + 1 ;
CE[13] := a[13]*(a[1] + a[3] ) + 1 ;
CE[14] := a[14]*(a[1]^2 + 2*a[2]) + 1 ;
CE[15] := a[15]*(a[2]*a[3] - a[4] + 1) + 1 ;
CE[16] := a[16]*a[4] + 1 ;

Za<[a]> := PolynomialRing(Integers(),16) ;
E := [] ;

for i in [1..16] do;
E[i] := Evaluate(CE[i], a ) ; 
end for ; 
 
CC := ComplexField(2000) ;
ZZa<[a]> := PolynomialRing(CC,16) ;
EE := [] ;
for i in [1..16] do ;
EE[i] := ZZa ! E[i] ;
end for ;

J := JacobianMatrix(EE) ;
im := CC.1 ;
e := Exp(1) ;
e := CC ! e ;



// the first approximation we use 

a1 := [ 2.89142578825226740000000000000 - 2.13725528220314360000000000000*CC.1, 
-0.435420544682335700000000000000 + 5.83909001972394000000000000000*CC.1, 
-4.57267582688547900000000000000 - 2.13725528220313600000000000000*CC.1, 
-4.27451056440629800000000000000 - 15.9526906038541640000000000000*CC.1, 
-28.7959743204547600000000000000 + 1.17995967957101000000000000000*CC.1, 
193.942454280126200000000000000 - 33.6619193888560800000000000000*CC.1, 
-215.824592091662480000000000000 + 219.824592091663330000000000000*CC.1, 
44.4544023365204600000000000000 - 195.762494600555360000000000000*CC.1, 
4.58005196605811050267170890417*(e^(-5)) + 4.30642082325788892079807918890*(e^(-5))*CC.1, 
-0.223652563063139160000000000000 - 0.165317271405363040000000000000*CC.1, 
0.179480926275014550000000000000 - 0.0838888590091154700000000000000*CC.1, 
-0.133974596215561500000000000000 - 1.57838174883301126209732825833*(e^(-7))*CC.1, 
0.0796874903416969500000000000000 - 0.202602237317134230000000000000*CC.1, 
-0.324623527238590950000000000000 - 0.0756936685597712900000000000000*CC.1, 
-0.0406074762465326240000000000000 - 0.0201893585838929070000000000000*CC.1, 
0.0156713375949008200000000000000 - 0.0584862281267338950000000000000*CC.1 ];

a1 := [ CC ! b : b in a1 ];
N := [a1] ;
J := JacobianMatrix(EE) ;

load "general.m";

// we perform 700 steps of Newton-Raphson 

for i in [2..700] do ; 
a := N[i-1] ;
a := [ CC ! b : b in a] ;
N[i] := nr(a,16) ;
end for ;

pt := N[700] ;


// we search for the minimal polynomial of  a[1] and relations for a[2], a[3] and a[4] in terms of a[1] 

m1 :=  minpoly(pt[1], 1000, 8);
r1 := re([pt[1],pt[2]], 700, 8) ;
r2 := re([pt[1], pt[3]], 700, 8) ;
r3 := re([pt[1], pt[4]], 700, 8) ;

// we verify that roots of the above equation are points on our scheme 

SK := SplittingField(m1) ;
R := Roots(m1,SK);
R1 := [ a[1] : a in R ];
R2 := [ [a] cat [Evaluate(r1,a),Evaluate(r2, a), Evaluate(r3,a)] : a in R1];

u1 := [ [-1/r[3], -1/(r[1] - r[3]), -1/(r[1] + r[3]), -1/( r[1]^2 + 2*r[2]), -1/(r[2]*r[3] - r[4] +1), -1/r[4] ] : r in R2 ]; 
r :=  [ re([pt[1],pt[i]],1000,8) : i in [5..10] ];
PT := [ [ Evaluate(s,a) : s in r ] : a in R1 ] ;
PT := [ R2[i] cat  PT[i] cat u1[i] : i in [1..#R] ] ;
EPT := [ [ Evaluate(e,a) : e in E ] : a in PT] ;
NET := [ a : a in EPT | a ne [ 0 : i in [1..16] ]] ;

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


// the second approximation we use is the following:

a1 := [0.43542054468233904 - 0.43542054468233904*im,  0.43542054468233904 + 0.43542054468233904*im, 0.24582949395087422 - 0.10703921215507343*im,   1.1895910507314649 + 0.31874996136678674*im,  -0.7499768091723923 + 0.16510367715278276*im, -0.19481568198918905 + 0.4953110314583483*im,   0.660414708611131 - 2.01159530831358*(e^-29)*im,      -0.9447924911615814 - 0.3302073543055655*im, 0.003191231510410138 - 0.007109196586347725*im, -1.1483151314432691 - 1.1483151314432691*im,  -3.4195450098619755 - 1.4889401507598758*im,    -1.3186276411015725 - 2.283930070652622*im,   -0.8983151314432691 - 0.7153024295510498*im, -0.8707603147552896 + 0.49161338421677636*im,   0.5284431270981033 - 3.7976676574353982*im,   -0.7843138205507862 + 0.21015625482915165*im];


a1 := [ CC ! b : b in a1 ];
N := [a1] ;
J := JacobianMatrix(EE) ;

load "general.m";

// we perform 700 steps of Newton-Raphson 

for i in [2..700] do ; 
a := N[i-1] ;
a := [ CC ! b : b in a] ;
N[i] := nr(a,16) ;
end for ;

pt := N[700] ;
// we search for the minimal polynomial of  a[1] and relations for a[2], a[3] and a[4] in terms of a[1] 

m1 :=  minpoly(pt[1], 1000, 8);
r1 := re([pt[1],pt[2]], 700, 8) ;
r2 := re([pt[1], pt[3]], 700, 8) ;
r3 := re([pt[1], pt[4]], 700, 8) ;



// we verify that roots of the above equation are points on our scheme 

SK := SplittingField(m1) ;
R := Roots(m1,SK);
R1 := [ a[1] : a in R ];
R2 := [ [a] cat [Evaluate(r1,a),Evaluate(r2, a), Evaluate(r3,a)] : a in R1];

u1 := [ [-1/r[3], -1/(r[1] - r[3]), -1/(r[1] + r[3]), -1/( r[1]^2 + 2*r[2]), -1/(r[2]*r[3] - r[4] +1), -1/r[4] ] : r in R2 ];


r :=  [ re([pt[1],pt[i]],700,8) : i in [5..6] ] ;
PT1 := [ [ Evaluate(s,a) : s in r ] : a in R1 ] ;
// the following expression for a[7] was found using factorisation over the number field defined by m1 
Zu<u> := PolynomialRing(Rationals());
s := 1/8*u^7 - 1/2*u^6 + 3/4*u^5 + u^4 - 5/2*u^3 + u;
PT2 := [ Evaluate(s,a) : a in R1];

rr2 := [ re([pt[1],pt[i]],1000,8) : i in [8..10] ] ;
PT3 := [ [ Evaluate(s,a) : s in rr2] : a in R1] ;

PT := [ PT1[i] cat [PT2[i]] cat PT3[i] : i in [1..#R] ];



PT := [ R2[i] cat  PT[i] cat u1[i] : i in [1..#R] ] ;
EPT := [ [ Evaluate(e,a) : e in E ] : a in PT] ;
NET := [ a : a in EPT | a ne [ 0 : i in [1..16] ]] ;

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

// the third orbit corresponds to the following approximation 

a1 := [ -0.1593749806833935 - 1.405204474634267*CC.1, 0.8406250193166056 + 
1.137255282203145*CC.1, -0.2885338913187149 - 0.5439947564300679*CC.1, 
-2.274510564406288 + 1.681250038633213*CC.1, 3.660254037844385 + 
3.729470089444420*CC.1, 1.586543195478669 + 2.316038140777203*CC.1, 
7.215390309173462 + 3.455161745314612*CC.1, -11.10320369387406 + 
4.868593693981836*CC.1, -0.007160349515405204 - 0.004277546691399913*CC.1, 
0.07968749034169680 - 0.7026022373171341*CC.1, 0.7609375289749105 - 
1.434653044886013*CC.1, -0.1703125096583028 - 1.135614939209353*CC.1, 
0.1119772180005272 - 0.4872998077660843*CC.1, 0.03580597972141455 + 
0.3637961893028698*CC.1, -0.1880642895110987 - 0.1270728045075998*CC.1, 
0.2843138205507863 + 0.2101562548291519*CC.1 ];

a1 := [ CC ! b : b in a1 ];
N := [a1] ;
J := JacobianMatrix(EE) ;

load "general.m";

// we perform 700 steps of Newton-Raphson 

for i in [2..700] do ; 
a := N[i-1] ;
a := [ CC ! b : b in a] ;
N[i] := nr(a,16) ;
end for ;

pt := N[700] ;
// we search for the minimal polynomial of  a[1] and relations for a[2], a[3] and a[4] in terms of a[1] 

m1 :=  minpoly(pt[1], 1000, 8);
r1 := re([pt[1],pt[2]], 700, 8) ;
r2 := re([pt[1], pt[3]], 700, 8) ;
r3 := re([pt[1], pt[4]], 700, 8) ;

// we verify that roots of the above equation are points on our scheme 

SK := SplittingField(m1) ;
R := Roots(m1,SK);
R1 := [ a[1] : a in R ];
R2 := [ [a] cat [Evaluate(r1,a),Evaluate(r2, a), Evaluate(r3,a)] : a in R1];

u1 := [ [-1/r[3], -1/(r[1] - r[3]), -1/(r[1] + r[3]), -1/( r[1]^2 + 2*r[2]), -1/(r[2]*r[3] - r[4] +1), -1/r[4] ] : r in R2 ];
r :=  [ re([pt[1],pt[i]],700,8) : i in [5..10] ];
PT := [ [ Evaluate(s,a) : s in r ] : a in R1 ] ;
PT := [ R2[i] cat  PT[i] cat u1[i] : i in [1..#R] ] ;
EPT := [ [ Evaluate(e,a) : e in E ] : a in PT] ;
NET := [ a : a in EPT | a ne [ 0 : i in [1..16] ]] ;

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


// the fourth orbit corresponds to the following approximation 

a1 := [ 1.456005243569932 - 1.159374980683393*CC.1, 0.7239544360010542 - 
0.1085742117477293*CC.1, 0.4560052435699316 + 0.5726758268854839*CC.1, 
3.464101615137754 + 1.362500077266427*CC.1, 3.583629154612597 + 
6.687811410608381*CC.1, 0.9836761543992252 - 4.280738613639845*CC.1, 
-1.569219381653055 + 5.984515691253272*CC.1, -0.6154373089601026 + 
7.976292178621589*CC.1, 0.0006901643118313084 - 0.003151744128679971*CC.1, 
-0.4203125096583033 - 0.3346827285946372*CC.1, -0.8509173687604029 + 
1.068627641101572*CC.1, -0.2500000000000000 - 0.4330127018922193*CC.1, 
-0.4780026217849658 - 0.1466747884494773*CC.1, -0.1245328840594731 - 
0.2012308878447313*CC.1, 0.3918523666380888 - 0.1886483808860123*CC.1, 
-0.2500000000000000 + 0.09832997329758218*CC.1 ];

a1 := [ CC ! b : b in a1 ];
N := [a1] ;
J := JacobianMatrix(EE) ;

load "general.m";

// we perform 700 steps of Newton-Raphson 

for i in [2..700] do ; 
a := N[i-1] ;
a := [ CC ! b : b in a] ;
N[i] := nr(a,16) ;
end for ;

pt := N[700] ;
// we search for the minimal polynomial of  a[1] and relations for a[2], a[3] and a[4] in terms of a[1] 

m1 :=  minpoly(pt[1], 1000, 8);
r1 := re([pt[1],pt[2]], 700, 8) ;
r2 := re([pt[1], pt[3]], 700, 8) ;
r3 := re([pt[1], pt[4]], 700, 8) ;

// we verify that roots of the above equation are points on our scheme 

SK := SplittingField(m1) ;
R := Roots(m1,SK);
R1 := [ a[1] : a in R ];
R2 := [ [a] cat [Evaluate(r1,a),Evaluate(r2, a), Evaluate(r3,a)] : a in R1];

u1 := [ [-1/r[3], -1/(r[1] - r[3]), -1/(r[1] + r[3]), -1/( r[1]^2 + 2*r[2]), -1/(r[2]*r[3] - r[4] +1), -1/r[4] ] : r in R2 ];
r :=  [ re([pt[1],pt[i]],700,8) : i in [5..10] ];
PT := [ [ Evaluate(s,a) : s in r ] : a in R1 ] ;
PT := [ R2[i] cat  PT[i] cat u1[i] : i in [1..#R] ] ;
EPT := [ [ Evaluate(e,a) : e in E ] : a in PT] ;
NET := [ a : a in EPT | a ne [ 0 : i in [1..16] ]] ;

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

// The last orbit of quadritangents used in our computation corresponds to the following approximation 

a1 := [ 0.384619775746675500000000000000 - 0.883329416684447600000000000000*im, 
1.79533990382431100000000000000 - 0.760731878024292400000000000000*im, 
-0.615380224253325300000000000000 + 0.848721390884429700000000000000*im, 
1.91459120526431860000000000000 - 0.217148423495457770000000000000*im, 
1.92150409572917600000000000000 + 2.08219076003366640000000000000*im, 
3.49518187665176330000000000000 - 3.90407850834409600000000000000*im, 
2.50773216421431180000000000000 - 1.60354614656269280000000000000*im, 
0.233943354409011650000000000000 + 0.216588024623247220000000000000*im, 
0.0348830235058634014502015755893 - 0.0508634200295954763349694479908*im, 
-0.414370218936334240000000000000 - 0.951655184848104500000000000000*im, 
0.559935001408853700000000000000 + 0.772252331925710900000000000000*im, 
-0.250000000000000060000000000000 - 0.433012701892219100000000000000*im, 
4.23817265096355200000000000000 - 0.635614939209352100000000000000*im, 
-0.217589265508115600000000000000 - 0.161882716669414500000000000000*im, 
0.203007562403547780000000000000 + 0.326440232047642550000000000000*im, 
-0.515671337594900600000000000000 - 0.0584862281267335350000000000000*im];

a1 := [ CC ! b : b in a1 ];
N := [a1] ;
J := JacobianMatrix(EE) ;

load "general.m";

// we perform 700 steps of Newton-Raphson 

for i in [2..700] do ; 
a := N[i-1] ;
a := [ CC ! b : b in a] ;
N[i] := nr(a,16) ;
end for ;


pt := N[700] ;
// we search for the minimal polynomial of  a[1] and relations for a[2], a[3] and a[4] in terms of a[1] 

m1 :=  minpoly(pt[1], 1000, 8);
r1 := re([pt[1],pt[2]], 700, 8) ;
r2 := re([pt[1], pt[3]], 700, 8) ;
r3 := re([pt[1], pt[4]], 700, 8) ;


// we verify that roots of the above equation are points on our scheme 

SK := SplittingField(m1) ;
R := Roots(m1,SK);
R1 := [ a[1] : a in R ];
R2 := [ [a] cat [Evaluate(r1,a),Evaluate(r2, a), Evaluate(r3,a)] : a in R1];

u1 := [ [-1/r[3], -1/(r[1] - r[3]), -1/(r[1] + r[3]), -1/( r[1]^2 + 2*r[2]), -1/(r[2]*r[3] - r[4] +1), -1/r[4] ] : r in R2 ];
r :=  [ re([pt[1],pt[i]],700,8) : i in [5..10] ];
PT := [ [ Evaluate(s,a) : s in r ] : a in R1 ] ;
PT := [ R2[i] cat  PT[i] cat u1[i] : i in [1..#R] ] ;
EPT := [ [ Evaluate(e,a) : e in E ] : a in PT] ;
NET := [ a : a in EPT | a ne [ 0 : i in [1..16] ]] ;

assert NET eq [] ;
print "One orbit of 4-tangent planes of the form -x[1] + a[1]x[2] + 
a[2]x[3] + a[3]x[4] + a[4]x[5] is as follows:";

print " The minimal polynomial of a[1] is:" ;
m1;

print " and relations for a[2], a[3], a[4] resp. in terms of a[1] are:";
r1;
r2;
r3;

O5 := [ m1, r1, r2, r3];

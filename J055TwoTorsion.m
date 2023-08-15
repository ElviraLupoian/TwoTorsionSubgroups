// This is the MAGMA code used to compute the 2-torsion subgroup of the J0(55) over the degree 12 number field stated in the paper


Z<U> := PolynomialRing(Integers());
f := 63607*U^12 + 868505*U^11 + 2747782*U^10 - 2704958*U^9 + 622467*U^8 + 245506*U^7 + 52935*U^6 - 18718*U^5 + 7263*U^4 + 4958*U^3 + 1212*U^2 + 65*U + 1;
Nf := NumberField(f) ;


K := OptimizedRepresentation(Nf) ;
OK := Integers(K) ;
PP := PrimesUpTo(50,K) ;
P := PP[4] ;

RR<[x]> := PolynomialRing(Rationals(),5) ;
f1 := x[1]*x[3] -x[2]^2 + x[2]*x[4] -x[2]*x[5] -x[3]^2 + 3*x[3]*x[4] + 
x[3]*x[5]-2*x[4]^2 -4*x[5]^2 ;
f2 := x[1]*x[4] -x[2]*x[3] +2*x[2]*x[4] -2*x[2]*x[5] -2*x[3]^2 + 4*x[3]*x[4] 
+5*x[3]*x[5] -2*x[4]^2 -4*x[4]*x[5] -3*x[5]^2 ;
f3 := x[1]*x[5] -2*x[2]*x[5] -x[3]^2 +2*x[3]*x[4] + x[3]*x[5] -x[4]^2 ;
P4 := ProjectiveSpace(RR) ;
X := Curve(P4, [f1,f2,f3] ) ;


F47,pi := ResidueClassField(P) ;
X47 := ChangeRing(X, F47) ;
Cl47, phi, psi := ClassGroup(X47) ;
Z := FreeAbelianGroup(1) ;
degr := hom< Cl47 -> Z | [ Degree(phi(g)) : g in OrderedGenerators(Cl47) ] > ;
J47 := Kernel(degr);

// above is the set of points over F_47 of J

// we define our orbits of quadritangents 
R := Roots(f,K) ;  

Z<u> := PolynomialRing(K);
S1 := 1/116604814277367729853897371860997284375*(-23716386961281248375413083410\
535652124941*u^11 - 322441000863133079344828351804700092922946*u^10 - 
1006800447077291694490737912187936888924827*u^9 + 
1051954185741553449977224514956563717422222*u^8 - 
342498183596167580885170604516263947246019*u^7 - 
21552206836977920785690319483445399654507*u^6 - 
30485121612436137224691102171925283899067*u^5 + 
5207593884792485252148797018362058349737*u^4 - 
4756952674783763423221596806109931364927*u^3 -
 1016896424158130618622042253960671076711*u^2 - 
403093366682588393265209109596980599882*u + 
10678424277043749924876401160721454368);
S2 := 1/116604814277367729853897371860997284375*(-40709882798252271615671755180\
242653635489*u^11 - 561954769107893605624695219552975204671784*u^10 -  
1841339426517125220287499996408579531349683*u^9 + 
1474472875003257177639417950031587884784613*u^8 - 
121578897624076939561791232414966906066901*u^7 - 
248460582705855998657814270694845599980128*u^6 - 
42852589059183734863644273672104098935543*u^5 + 
7027290777565567354434871044406544298523*u^4 - 
3208548051200430712501289477244462464733*u^3 -
    4499554506820319627791342778895783365069*u^2 - 
1151774282830835009372954827072887837003*u - 
23026966226537449715903277566750139303);
S3 := 1/23320962855473545970779474372199456875*(6873991845110251896204725474482\
3217879754*u^11 + 932620077111546200919511530409190847671959*u^10 + 
2888579126432963149439181977471551862302808*u^9 - 
3173103576562428862198101093730255975006553*u^8 + 
    952763992879910830668680212088814158948331*u^7 + 
184232668488364426267088779799237874639493*u^6 + 
    37413053808909232030679933383018974146803*u^5 - 
22566976960190982309236680802672758999803*u^4 + 
10516701627534286288292865217890722266918*u^3
    + 4598609562394840488297318272662632643084*u^2 + 
855607141243614589571098805261264899518*u - 
25283332426994990717802161322094001062);
coeff1 := [ [R[i,1], Evaluate(S1, R[i,1]), Evaluate(S2,R[i,1]),-1,  Evaluate(S3,
R[i,1]) ] : i in [1..6] ];
R1 := R ;

f6 := 144057127*U^12 + 27502925*U^11 + 98022670*U^10 + 55316502*U^9 + 
5082695*U^8 + 306018*U^7 + 352531*U^6 - 75086*U^5 - 19013*U^4 + 2234*U^3 + 180*U^2 - 23*U + 1;
R := Roots(f6,K) ;
Z<u> := PolynomialRing(K);
S1 := 1/5319915865418302576715947994898416108819611888843*(-4756712737048165633\
5821211712562509108785823475718599340*u^11 - 
    4323070432008161576376898646208322714675319745984659241*u^10 - 
    32603641813335708511197118593314062202961475975653547119*u^9 - 
    15079654281880675224966830906094965718904700575159253850*u^8 - 
    626713696085727658459663304605159820221581572798883153*u^7 - 
    261243838955520648703956935107953430888018561700643341*u^6 - 
    101497623650378851454658228955834691344870258064944283*u^5 + 
    31536613428598963794001394673206711350314332628420595*u^4 + 
    1544385768325498707275776679060165635746011801485355*u^3 - 
    486888719912119609387888101841651842577902839884951*u^2 + 
    11168686935604908625963408156935851283136207335393*u - 
    2350162398453525672300414629158870769619423241951);
S2 := 1/5319915865418302576715947994898416108819611888843*(56066411217915833770\
002228948021751021711032886178865844*u^11 +  
2196409684859593221020442978167230035496057397820055386*u^10 +  
38094019816882513763280782423474600100280487310544216973*u^9 + 
15783596246579429073488904845893037286536198547019708974*u^8 - 
228081579797602085556690048651449663402689075475596257*u^7 + 
    248420803515931574949740736332938408873359448374641855*u^6 + 
    104791055542528130825365214413560001400311074204547169*u^5 - 
    43992555761481046880738613044998250443532516088075708*u^4 - 
    63608951887866718724494294762241979088860107679277*u^3 + 
    719417092575540633538098459014647977809787974632212*u^2 - 
    54063648941202881167941087774049199551368942651638*u + 
    7825433308483606644305239944727695181518369986884);
S3 := 1/5319915865418302576715947994898416108819611888843*(23676709321084869484\
525441437108183403041091361936447497*u^11 - 
20292626808214509199368399658785226402393662647693860765*u^10 + 
14191135396498134745004697915476997824103534660785779742*u^9 - 
7865378322631514839488089503049699445731491354302575768*u^8 - 
6805815658498029699549085697746901745889885113509186650*u^7 - 
    160236371915362333847073153713086504272011640402569303*u^6 - 
    71647096993487269594621830836266016601225124646267452*u^5 - 
    67681536791625696356733832192517684492476917290826354*u^4 + 
    13831829866788193841205126104664244395948218248731174*u^3 + 
    935966892562957672450862462100602393651386851056728*u^2 - 
    291710536101541940039070198967174725583469259362154*u + 
    18873055854733158679886230860731243116156729292642);

coeff6 := [ [R[i,1], Evaluate(S1, R[i,1]), Evaluate(S2,R[i,1]),-1,  Evaluate(S3,
R[i,1]) ] : i in [1..6] ];
R2 := R ;

f11 := 6737*U^6 - 3838*U^5 - 3249*U^4 + 1604*U^3 + 1527*U^2 - 902*U + 121;
R := Roots(f11,K) ;
Z<u> := PolynomialRing(K);
S1 := (1/98850400)*(19127939669*u^5 - 6837632467*u^4 - 10846801480*u^3 + 
2393348508*u^2 +  4886923427*u - 1361707457);
S2 := u + 1 ;
S3 := (1/49425200)*(43352851006*u^5 - 10993587083*u^4 - 24615807720*u^3 + 
2580823642*u^2 +  10609943598*u - 2597943843);
coeff11 := [ [R[i,1], Evaluate(S1, R[i,1]), Evaluate(S2,R[i,1]),-1,  Evaluate(S3, R[i,1]) ] : i in [1..6] ];
R3 := R ;

coeff := coeff1 cat coeff6 cat coeff11; 

// we clear denominators and scale to ensure all planes are defined over the ring of integer OK


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

FX47 := FunctionField(X47) ;
CQT := [ [MonomialCoefficient(a,X[1]),  MonomialCoefficient(a,X[2]) , MonomialCoefficient(a,X[3]) , MonomialCoefficient(a,X[4]), 
MonomialCoefficient(a,X[5])] : a in QT ] ;
F47QT := [ [ pi(s) : s in a ] : a in CQT ] ;

ZU<[U]> := PolynomialRing(F47, 5) ;

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
hh := hom< ZN -> J47 | [ a : a in H ] > ;
ihh := Image(hh) ;

print "The 2-torsion subgroup generated by our orbits of quadritangents is";
ihh;

print "Since the set of F_47 rational points of J is:" ;
J47 ;
print "and reduction modulo 47 is injective when restristed to the K-rational torsion points:";
print "our subgroup is the J(K)[2], where K is";
K;
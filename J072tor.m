

Zx<[x]> := PolynomialRing(Integers(), 5) ; 

f1 := x[1]*x[3] -x[2]^2 + x[4]^2 ;  
f2 :=  x[1]*x[4] - x[3]^2 ;
f3 := x[1]*x[5] - x[3]*x[4] -2*x[5]^2 ;


C2 := Curve(ProjectiveSpace(Zx), [f1, f2, f3]);


// the rational cusps 

T := [[1,0,0,0,0], [1,0,-1, 1,1],  [0,1,0,1,0], [0,-1, 0, 1,0], [-2, 0, 2,-2, 1], [4,3, 2, 1, 1], [2, 0,0,0,1],[4, -3, 2, 1,1] ];

X := ChangeRing(C2, GF(5));
Cl, phi, psi := ClassGroup(X);
Z := FreeAbelianGroup(1) ; 
degr := hom<Cl -> Z | [Degree(phi(g)) : g in OrderedGenerators( Cl)]>;
J := Kernel(degr);

pts := [ X ! a : a in T];
pts := [ Place(a) : a in pts];

D1 := -5*pts[7] + 2*pts[6] + pts[8] + pts[3] + pts[1] ;
D2 := 3*pts[5] + 3*pts[3] + 6*pts[1] -12*pts[7];
D3 := pts[2] - pts[1];

ZN := FreeAbelianGroup(3) ;
hh := hom< ZN -> J | [ psi(D1),psi(D2), psi(D3)] >;
H := Image(hh);  

assert #H eq 4^3 ;   // we verify that the above generate the rational 4 torsion subgroup 


D := -psi(D1) + 2*psi(D2) + psi(D3);
P := 12*J.4 ;  // this is an 8 torsion point mod 5

// if we assume J0(72)(Q) =  Z/2 x Z/4 x Z/12 x Z/24 then P is necessarily the reduction of a rational 8 torsion point (by injectivity of torsion 
// we find its cuspdial representation 

assert  2*P eq D ; // thus if a rational 8 torsion point exists, then there exists such a point such that twice it = D 

// we reduce mod 7 
T := [[1,0,0,0,0], [1,0,-1, 1,1],  [0,1,0,1,0], [0,-1, 0, 1,0], [-2, 0, 2,-2, 1], [4,3, 2, 1, 1], [2, 0,0,0,1],[4, -3, 2, 1,1] ];
X := ChangeRing(C2, GF(7));
Cl, phi, psi := ClassGroup(X);
Z := FreeAbelianGroup(1) ; 
degr := hom<Cl -> Z | [Degree(phi(g)) : g in OrderedGenerators( Cl)]>;
J := Kernel(degr);


pts := [ X ! a : a in T];
pts := [ Place(a) : a in pts];

D1 := -5*pts[7] + 2*pts[6] + pts[8] + pts[3] + pts[1] ;
D2 := 3*pts[5] + 3*pts[3] + 6*pts[1] -12*pts[7];
D3 := pts[2] - pts[1];

ZN := FreeAbelianGroup(3) ;
hh := hom< ZN -> J | [ psi(D1),psi(D2), psi(D3)] >;
H := Image(hh);   // image of the rational 4 torsion  mod 7
assert #H eq 4^3 ;

ZN := FreeAbelianGroup(2);
h2 := hom< ZN -> J | [ 6*J.5, 6*J.6] >;
H2 := Image(h2) ;

H3 := H2 meet H ; // the  group of 4 torsion points which are cuspidal and twice a rational 8 torsion point 


D := -psi(D1) + 2*psi(D2) + psi(D3); // this is twice times a rational 8 torsion point (if one exits)

TT := [ a : a in H3 | a eq D];
TT ;



/*
  Parameterized exact test for the candidate

      Re(z1^2*conj(z2)) + eps*|z1|^2*Im(z1^2*conj(z2))

  in any fixed dimension d >= 4.  Change d below and run Magma.

  A successful run prints that the saturated ideal is the unit ideal.  That
  output is a rigorous certificate for the selected dimension.  The system is
  quadratic/cubic but grows with d; no claim is made here that every selected
  dimension will finish quickly.
*/

d := 5;
assert d ge 4;
Nverts := d + 1;

/* Number of variables after eliminating all zero sums and fixing Im(a0)=0:
     4*Nverts-5 coordinate variables,
     9 exactness parameters,
     6 auxiliary third moments,
     1 inverse variable. */
nvars := 4*Nverts + 11;
Q := Rationals();
P := PolynomialRing(Q, nvars, "grevlex");

names := [];
names cat:= [ "x" cat IntegerToString(i) : i in [0..Nverts-2] ];
names cat:= [ "y" cat IntegerToString(i) : i in [1..Nverts-2] ];
names cat:= [ "p" cat IntegerToString(i) : i in [0..Nverts-2] ];
names cat:= [ "q" cat IntegerToString(i) : i in [0..Nverts-2] ];
names cat:= [ "sigR", "sigI", "m", "nuR", "nuI",
              "ellR", "ellI", "lamR", "lamI",
              "vvR", "vvI", "kkR", "kkI", "rhoR", "rhoI",
              "zinv" ];
assert #names eq nvars;
AssignNames(~P, names);
z := [ P.i : i in [1..nvars] ];

pos := 1;
xshort := z[pos..pos+Nverts-2]; pos +:= Nverts-1;
yshort := z[pos..pos+Nverts-3]; pos +:= Nverts-2;
pshort := z[pos..pos+Nverts-2]; pos +:= Nverts-1;
qshort := z[pos..pos+Nverts-2]; pos +:= Nverts-1;

sig := [z[pos],z[pos+1]]; pos +:= 2;
m := z[pos]; pos +:= 1;
nu := [z[pos],z[pos+1]]; pos +:= 2;
ell := [z[pos],z[pos+1]]; pos +:= 2;
lam := [z[pos],z[pos+1]]; pos +:= 2;
vv := [z[pos],z[pos+1]]; pos +:= 2;
kk := [z[pos],z[pos+1]]; pos +:= 2;
rho := [z[pos],z[pos+1]]; pos +:= 2;
zinv := z[pos]; pos +:= 1;
assert pos eq nvars+1;

x := xshort cat [ -&+xshort ];
y := [ P!0 ] cat yshort cat [ -&+yshort ];
p := pshort cat [ -&+pshort ];
q := qshort cat [ -&+qshort ];

function CAdd(u, v)
  return [u[1]+v[1],u[2]+v[2]];
end function;

function CSub(u, v)
  return [u[1]-v[1],u[2]-v[2]];
end function;

function CScale(r, u)
  return [r*u[1],r*u[2]];
end function;

function CMul(u, v)
  return [u[1]*v[1]-u[2]*v[2],u[1]*v[2]+u[2]*v[1]];
end function;

function CConj(u)
  return [u[1],-u[2]];
end function;

function CAbs2(u)
  return u[1]^2+u[2]^2;
end function;

function CSum(v)
  return [&+[u[1]:u in v],&+[u[2]:u in v]];
end function;

function CProd(v)
  ans := [P!1,P!0];
  for u in v do
    ans := CMul(ans,u);
  end for;
  return ans;
end function;

a := [[x[i],y[i]]:i in [1..Nverts]];
b := [[p[i],q[i]]:i in [1..Nverts]];

eqs := [];

/* Real Stiefel-frame equations for a,b in 1^perp. */
eqs cat:= CSum([CMul(a[i],a[i]):i in [1..Nverts]]);
eqs cat:= CSum([CMul(b[i],b[i]):i in [1..Nverts]]);
eqs cat:= CSum([CMul(a[i],b[i]):i in [1..Nverts]]);
eqs cat:= CSum([CMul(a[i],CConj(b[i])):i in [1..Nverts]]);
Append(~eqs,&+[CAbs2(a[i]):i in [1..Nverts]]-2);
Append(~eqs,&+[CAbs2(b[i]):i in [1..Nverts]]-2);
assert #eqs eq 10;

/* Universal exactness identities.  They are equivalent to exact Q-labels:

   a^2 = sig*a + (m/2)b + (nu/2)conj(a) + ell*conj(b),
   a*conj(b) = lam*a + (conj(nu)/4)b
                         + (m/2)conj(a) + (sig/2)conj(b).
*/
for i in [1..Nverts] do
  rhs1 := CSum([
    CMul(sig,a[i]),
    CScale(m/2,b[i]),
    CScale(1/2,CMul(nu,CConj(a[i]))),
    CMul(ell,CConj(b[i]))
  ]);
  eqs cat:= CSub(CMul(a[i],a[i]),rhs1);

  rhs2 := CSum([
    CMul(lam,a[i]),
    CScale(1/4,CMul(CConj(nu),b[i])),
    CScale(m/2,CConj(a[i])),
    CScale(1/2,CMul(sig,CConj(b[i])))
  ]);
  eqs cat:= CSub(CMul(a[i],CConj(b[i])),rhs2);
end for;

/* Auxiliary moments used to express the obstruction cubically. */
vvdef := CSum([CMul(a[i],CMul(b[i],b[i])):i in [1..Nverts]]);
kkdef := CSum([CProd([b[i],b[i],b[i]]):i in [1..Nverts]]);
rhodef := CSum([CMul(CMul(b[i],b[i]),CConj(b[i])):i in [1..Nverts]]);
eqs cat:= CSub(vv,vvdef) cat CSub(kk,kkdef) cat CSub(rho,rhodef);

/* If E denotes the following complex expression, then on the exactness
   locus Re(E)=8*sqrt(2)*Lambda(Y_H).  This identity follows by expanding the
   degree-five edge expression in wedges and reducing with the two identities
   above. */
mp := [m,P!0];
Eterms := [
  CScale(-3/2,CProd([mp,mp,mp])),
  CScale( 3/4,CProd([CConj(nu),mp,nu])),
  CScale(-12,CProd([CConj(rho),ell,lam])),
  CScale(-16,CProd([CConj(lam),lam,mp])),
  CScale(-2,CProd([CConj(nu),CConj(sig),ell])),
  CScale(-6,CProd([rho,lam,mp])),
  CScale(-8,CProd([CConj(sig),lam,nu])),
  CScale(-8,CProd([lam,sig,sig])),
  CScale(-1,CProd([CConj(vv),mp,vv])),
  CScale(12,CProd([CConj(kk),CConj(lam),ell])),
  CScale(14,CProd([CConj(ell),ell,mp])),
  CScale(16,CProd([ell,lam,lam])),
  CScale(4,CProd([CConj(lam),CConj(vv),nu])),
  CScale(5,CProd([CConj(ell),nu,sig])),
  CScale(6,CProd([CConj(lam),CConj(rho),mp])),
  CScale(6,CProd([CConj(vv),ell,sig]))
];
E := CSum(Eterms);

/* Gauge-chart saturation and the forbidden zero obstruction. */
Append(~eqs,zinv*x[1]-1);
Append(~eqs,E[1]);

I := ideal<P|eqs>;
SetGBGlobalModular(true);
print "dimension",d,"variables",nvars,"equations",#eqs;
time proper := IsProper(I);
assert not proper;
print "Certificate verified in dimension",d,
      ": the saturated ideal is the unit ideal.";


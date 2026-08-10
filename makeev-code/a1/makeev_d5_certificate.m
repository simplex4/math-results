/*
  Exact algebraic certificate for the proposed A_5 counterexample.

  If this program terminates with

      Certificate verified: the saturated ideal is the unit ideal.

  then the algebraic-transversality assertion in makeev_d5_note.tex is
  proved over C (and hence over R).  All coefficients are rational.

  Tested syntax target: Magma V2.29.  The expensive line is IsProper(I).
*/

Q := Rationals();
P := PolynomialRing(Q, 25, "grevlex");
AssignNames(~P,
  [ "x0", "x1", "x2", "x3", "x4",
    "y1", "y2", "y3", "y4",
    "s0", "s1", "s2", "s3", "s4",
    "t0", "t1", "t2", "t3", "t4",
    "c0", "c1", "c2", "c3", "c4", "zinv" ]);
z := [ P.i : i in [1..25] ];

/* Eliminate the zero-sum coordinate in each real vector.  The circle gauge
   is y0=0.  Thus a_i=x_i+i*y_i and b_i=s_i+i*t_i, for i=0,...,5. */
xshort := z[1..5];
yshort := z[6..9];
sshort := z[10..14];
tshort := z[15..19];
cshort := z[20..24];
zinv := z[25];

x := xshort cat [ -&+xshort ];
y := [ P!0 ] cat yshort cat [ -&+yshort ];
s := sshort cat [ -&+sshort ];
t := tshort cat [ -&+tshort ];
c := cshort cat [ -&+cshort ];

/* Complex numbers are represented as [real part, imaginary part]. */
function CAdd(u, v)
  return [ u[1] + v[1], u[2] + v[2] ];
end function;

function CSub(u, v)
  return [ u[1] - v[1], u[2] - v[2] ];
end function;

function CScale(r, u)
  return [ r*u[1], r*u[2] ];
end function;

function CMul(u, v)
  return [ u[1]*v[1] - u[2]*v[2],
           u[1]*v[2] + u[2]*v[1] ];
end function;

function CConj(u)
  return [ u[1], -u[2] ];
end function;

function CAbs2(u)
  return u[1]^2 + u[2]^2;
end function;

function CSum(v)
  return [ &+[ u[1] : u in v ], &+[ u[2] : u in v ] ];
end function;

function CPow(u, n)
  ans := [ Parent(u[1])!1, Parent(u[1])!0 ];
  for j in [1..n] do
    ans := CMul(ans, u);
  end for;
  return ans;
end function;

a := [ [x[i], y[i]] : i in [1..6] ];
b := [ [s[i], t[i]] : i in [1..6] ];

eqs := [];

/* Frame equations.  Together with the eliminated sums, these say that
   Re(a), Im(a), Re(b), Im(b) are an orthonormal real 4-frame in 1^perp. */
aa := CSum([ CMul(a[i], a[i]) : i in [1..6] ]);
bb := CSum([ CMul(b[i], b[i]) : i in [1..6] ]);
ab := CSum([ CMul(a[i], b[i]) : i in [1..6] ]);
aB := CSum([ CMul(a[i], CConj(b[i])) : i in [1..6] ]);
eqs cat:= aa cat bb cat ab cat aB;
Append(~eqs, &+[ CAbs2(a[i]) : i in [1..6] ] - 2);
Append(~eqs, &+[ CAbs2(b[i]) : i in [1..6] ] - 2);

/* c is the unit real vector completing that 4-frame inside 1^perp. */
Append(~eqs, &+[ c[i]^2 : i in [1..6] ] - 1);
Append(~eqs, &+[ c[i]*x[i] : i in [1..6] ]);
Append(~eqs, &+[ c[i]*y[i] : i in [1..6] ]);
Append(~eqs, &+[ c[i]*s[i] : i in [1..6] ]);
Append(~eqs, &+[ c[i]*t[i] : i in [1..6] ]);

/* Exactness of the Q-labels, in the nine-equation form derived in the note.
   m is real; 2 Sum(a*conj(b)^2)=conj(Sum(a^3)); the components of a^2 and
   a*conj(b) along c vanish; and 2 Sum(a*|b|^2)=Sum(|a|^2*a). */
m := CSum([ CMul(CMul(a[i], a[i]), CConj(b[i])) : i in [1..6] ]);
Append(~eqs, m[2]);

twou := CScale(2, CSum([
  CMul(a[i], CMul(CConj(b[i]), CConj(b[i]))) : i in [1..6]
]));
nbar := CConj(CSum([ CPow(a[i], 3) : i in [1..6] ]));
un := CSub(twou, nbar);
eqs cat:= un;

rc := CSum([ CScale(c[i], CMul(a[i], a[i])) : i in [1..6] ]);
wc := CSum([ CScale(c[i], CMul(a[i], CConj(b[i]))) : i in [1..6] ]);
eqs cat:= rc cat wc;

dleft := CScale(2, CSum([
  CScale(CAbs2(b[i]), a[i]) : i in [1..6]
]));
dright := CSum([
  CScale(CAbs2(a[i]), a[i]) : i in [1..6]
]);
eqs cat:= CSub(dleft, dright);

assert #eqs eq 24;

/* L is 4*sqrt(2) times Lambda_g(Y_H(g)); its zero set is rational.
   For da=a_i-a_j and db=b_i-b_j,
     Im(da^2*conj(db))
   is the unnormalized R-label. */
L := P!0;
for i in [1..6] do
  for j in [i+1..6] do
    da := CSub(a[i], a[j]);
    db := CSub(b[i], b[j]);
    xij := a[i][2]*a[j][1] - a[i][1]*a[j][2]
            + 2*(b[i][2]*b[j][1] - b[i][1]*b[j][2]);
    rlabel := CMul(CMul(da, da), CConj(db))[2];
    L +:= xij*CAbs2(da)*rlabel;
  end for;
end for;

/* Every circle orbit can, after a vertex permutation, be put in the chart
   y0=0, x0!=0.  The inverse variable performs the saturation by x0. */
Append(~eqs, zinv*x[1] - 1);
Append(~eqs, L);
assert #eqs eq 26;

I := ideal< P | eqs >;
SetGBGlobalModular(true);
time proper := IsProper(I);

assert not proper;
print "Certificate verified: the saturated ideal is the unit ideal.";


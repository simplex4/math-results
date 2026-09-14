# Prompt

Is there a cubic polynomial $f(x,y)$ with integer coefficients that is injective $\mathbb Z^2\to\mathbb Z$? (I know that an injective quartic exists, and no injective quadratics exist.)

# Output

**No—a Pell-equation argument rules out cubics.** I have not found a reference for the full argument below; the discussions I found leave the general, nonhomogeneous case unanswered. Thus, what follows is an independent derivation, rather than a quotation of an established result. ([math.stackexchange.com][1])

The key is to look for collisions at two integer points symmetric about a suitably chosen midpoint.

## 1. A lemma using Pell’s equation

**Lemma.** Let $q$ be a nondegenerate indefinite binary quadratic form with integer coefficients, let $\ell$ be a nonzero integer linear form, and let $c\in\mathbb Z\setminus\{0\}$. Then

$$
q(z)+\ell(z)+ct^2=0
$$

has a solution $z\in\mathbb Z^2$, $t\in\mathbb Z\setminus\{0\}$.

**Proof.** We first choose $u\in\mathbb Z^2$ such that, writing

$$
a=q(u),\qquad b=\ell(u),
$$

we have

$$
b\ne0,\qquad D=-ac>0,\qquad D\text{ is not a square}.
\tag{1}
$$

Here is why that choice is possible. Since $q$ is indefinite, choose an integer vector $r$ with $-cq(r)>0$ and $\ell(r)\ne0$, and an integer vector $s$ independent of $r$. The polynomial

$$
R(n)=-cq(nr+s)
$$

is quadratic, with positive leading coefficient and nonzero discriminant. It cannot be a square for every sufficiently large integer $n$: otherwise, writing $R(n)=m_n^2$ with $m_n>0$, the integer differences $m_{n+1}-m_n$ would converge to the square root of the leading coefficient. They would therefore eventually be constant, forcing $R$ itself to be the square of a linear polynomial, contrary to its nonzero discriminant. Thus $u=nr+s$ satisfies (1) for a suitable sufficiently large $n$.

Pell’s equation supplies integers $U>1$, $V>0$ with

$$
U^2-DV^2=1,
\quad\text{equivalently}\quad
U^2+acV^2=1.
$$

This is the standard existence theorem for Pell’s equation with positive nonsquare $D$. ([Lean Community][2])

Set

$$
k=-bcV^2,\qquad t=bUV,\qquad z=ku.
$$

Then $t\ne0$, and

$$
\begin{aligned}
q(z)+\ell(z)+ct^2
&=ak^2+bk+ct^2\\
&=b^2cV^2\bigl(acV^2-1+U^2\bigr)\\
&=0.
\end{aligned}
$$

This proves the lemma. $\square$

## 2. Cubics whose leading part is not a cube of a linear form

Let $H$ denote the homogeneous degree-three part of $f$. For a vector $v$, write

$$
\partial_v f=v\cdot\nabla f.
$$

The exact central-difference identity for a cubic is

$$
\boxed{\;
f(z+tv)-f(z-tv)
=2t\bigl(\partial_v f(z)+t^2H(v)\bigr).
\;}
\tag{2}
$$

Suppose first that $H$ is **not** a scalar multiple of the cube of a linear form.

We will find integer vectors $p,v$, with $v\ne0$, for which

$$
\partial_v f(p)=0,
$$

the quadratic part of $\partial_v f$ is nondegenerate and indefinite, and

$$
\nabla(\partial_v f)(p)\ne0.
\tag{3}
$$

### Choosing the point and direction

Let $M_H(w)$ denote the Hessian matrix of $H$ at $w$. The determinant of this Hessian is not identically zero. Indeed, if

$$
H(x,y)=a x^3+b x^2y+cxy^2+dy^3,
$$

then

$$
\det M_H(x,y)
=
4\bigl((3ac-b^2)x^2+(9ad-bc)xy+(3bd-c^2)y^2\bigr),
$$

and setting all three coefficients to zero forces $H$ to be a scalar cube of a rational linear form.

Consequently, choose $w\in\mathbb Z^2$ with $M_H(w)$ invertible. Put

$$
J(\alpha,\beta)=(\beta,-\alpha),
\qquad
v_0=J\nabla H(w).
$$

The vector $v_0$ is nonzero: Euler’s identity gives

$$
M_H(w)w=2\nabla H(w).
$$

Consider the binary quadratic form

$$
q_0(z)=\partial_{v_0}H(z).
$$

It satisfies

$$
q_0(w)=v_0\cdot\nabla H(w)=0,
$$

but

$$
\nabla q_0(w)=M_H(w)v_0\ne0.
$$

A real binary quadratic form with a nonzero zero at which its gradient is nonzero must be nondegenerate and indefinite. Thus $q_0$ has those properties.

Now take

$$
p_N=Nw,\qquad v_N=J\nabla f(p_N).
$$

These are integer vectors, and automatically

$$
\partial_{v_N}f(p_N)=0.
$$

As $N\to\infty$,

$$
\frac{v_N}{N^2}\longrightarrow v_0,
\qquad
\frac{M_f(p_N)}{N}\longrightarrow M_H(w).
$$

Therefore, for all sufficiently large integer $N$, the form

$$
q(z)=\partial_{v_N}H(z)
$$

is nondegenerate and indefinite, while $v_N\ne0$ and $M_f(p_N)$ is invertible.

Fix such an $N$, and abbreviate $p=p_N$, $v=v_N$. Then

$$
\nabla(\partial_v f)(p)=M_f(p)v\ne0,
$$

so all the conditions in (3) hold.

### Producing the collision

Since $\partial_v f(p)=0$, translating this quadratic polynomial gives

$$
\partial_v f(p+z)=q(z)+\ell(z),
$$

where

$$
q(z)=\partial_v H(z),
\qquad
\ell(z)=z^{\mathsf T}M_f(p)v.
$$

Here $q$ is nondegenerate and indefinite, and $\ell\ne0$.

Set $c=H(v)$. If $c=0$, equation (2) immediately gives

$$
f(p+v)=f(p-v).
$$

If $c\ne0$, the lemma provides integer $z$ and nonzero integer $t$ such that

$$
\partial_v f(p+z)+t^2H(v)=0.
$$

Applying (2) again gives

$$
\boxed{\;
f(p+z+tv)=f(p+z-tv).
\;}
$$

The two inputs are distinct because $t\ne0$ and $v\ne0$. This rules out injectivity whenever $H$ is not a scalar cube.

## 3. The remaining case: the leading part is a scalar cube

Suppose $H$ is a scalar multiple of the cube of a rational linear form. A unimodular integer change of variables puts $H$ in the form $a x^3$. Such a change preserves injectivity on $\mathbb Z^2$.

The whole polynomial then has the form

$$
f(x,y)=Ay^2+B(x)y+P(x),
\tag{4}
$$

where $A\in\mathbb Z$ and $B,P\in\mathbb Z[x]$. We can give elementary collisions for every polynomial of this form.

Write $b=B(0)$. If $b=0$, then

$$
f(0,1)=f(0,-1).
$$

If $b\ne0$ but $A=0$, then

$$
f\!\left(0,\frac{P(b)-P(0)}b\right)=f(b,0),
$$

and the displayed quotient is an integer.

It remains to handle $Ab\ne0$. Completing the square in (4), define

$$
R(x)=4AP(x)-B(x)^2,
$$

so that

$$
4Af(x,y)=(2Ay+B(x))^2+R(x).
$$

Set

$$
r=8Ab,\qquad
m=\frac{R(r)-R(0)}{8Ab},
\qquad
n=\frac{b-B(r)}{2A}-m.
$$

Both $m,n$ are integers: an integer polynomial satisfies
$S(r)\equiv S(0)\pmod r$, and $2A\mid r$.

Now

$$
\begin{aligned}
4A\bigl(f(0,m)-f(r,n)\bigr)
&=(b+2Am)^2-(b-2Am)^2+R(0)-R(r)\\
&=8Abm+R(0)-R(r)\\
&=0.
\end{aligned}
$$

Thus

$$
f(0,m)=f(r,n),
$$

with distinct inputs since $r\ne0$.

This also excludes the scalar-cube case.

---

Hence every integer-coefficient cubic in two variables has an integer collision. Together with the quartic construction and the quadratic obstruction you mentioned, **the minimum possible degree is $4$**.

[1]: https://math.stackexchange.com/questions/1073961/is-there-an-injective-cubic-polynomial-mathbb-z2-rightarrow-mathbb-z "number theory - Is There An Injective Cubic Polynomial $\mathbb Z^2 \rightarrow \mathbb Z$? - Mathematics Stack Exchange"
[2]: https://leanprover-community.github.io/mathlib4_docs/Mathlib/NumberTheory/Pell.html?utm_source=chatgpt.com "Mathlib.NumberTheory.Pell"

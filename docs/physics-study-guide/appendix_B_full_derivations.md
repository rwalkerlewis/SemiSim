# Appendix B — Full derivations

This appendix collects the long algebraic derivations that the main
chapters reference but do not unpack in detail. Skip on first read;
return when you need to verify a step.

## B.1 — Effective densities of states

Goal: derive

$$
N_c = 2\left(\frac{m_n^* kT}{2\pi\hbar^2}\right)^{3/2}.
\tag{B.1}
$$

Start from (3.2):

$$
n = \int_{E_c}^{\infty} g_c(E)\,f_\mathrm{FD}(E; E_F, T)\,dE.
$$

The 3D parabolic-band density of states is (2.4):

$$
g_c(E) = \frac{1}{2\pi^2}\left(\frac{2m_n^*}{\hbar^2}\right)^{3/2}\sqrt{E - E_c}.
$$

Multiply by 2 for spin and substitute. In the Boltzmann limit
$f_\mathrm{FD} \approx \exp(-(E - E_F)/kT)$:

$$
n = \frac{2}{2\pi^2}\left(\frac{2m_n^*}{\hbar^2}\right)^{3/2}\,e^{(E_F - E_c)/kT}
   \int_{E_c}^\infty \sqrt{E - E_c}\,e^{-(E - E_c)/kT}\,dE.
$$

Substitute $\eta = (E - E_c)/kT$, $d\eta = dE/kT$:

$$
n = \frac{1}{\pi^2}\left(\frac{2m_n^*}{\hbar^2}\right)^{3/2}\,(kT)^{3/2}\,e^{(E_F - E_c)/kT}
   \int_0^\infty \sqrt\eta\,e^{-\eta}\,d\eta.
$$

The integral is $\Gamma(3/2) = \sqrt\pi/2$. Combining:

$$
n = \frac{\sqrt\pi}{2\pi^2}\left(\frac{2m_n^* kT}{\hbar^2}\right)^{3/2}\,e^{(E_F-E_c)/kT}
   = \frac{1}{2\pi^{3/2}}\left(\frac{2m_n^* kT}{\hbar^2}\right)^{3/2}\,e^{(E_F-E_c)/kT}.
$$

Identifying $n = N_c\exp(-(E_c - E_F)/kT)$:

$$
N_c = \frac{1}{2\pi^{3/2}}\left(\frac{2m_n^* kT}{\hbar^2}\right)^{3/2}
    = 2\left(\frac{m_n^* kT}{2\pi\hbar^2}\right)^{3/2},
$$

using $(2/\hbar^2)^{3/2}\cdot 1/(2\pi^{3/2}) = 1/(\pi^{3/2}\hbar^3)\cdot 2^{1/2}
= 2/(2\pi\hbar^2/m^*kT)^{3/2}$ after rearrangement. (B.1) ✓

The hole expression $N_v = 2(m_p^*kT/(2\pi\hbar^2))^{3/2}$ follows by
the symmetric integral over the valence band.

## B.2 — Full SRH derivation

Goal: derive (6.5):

$$
R_\mathrm{SRH}(n, p) = \frac{np - n_i^2}{\tau_p(n + n_1) + \tau_n(p + p_1)}.
$$

Steady-state trap occupancy (6.2):

$$
f_t = \frac{c_n n + e_p}{c_n n + e_p + c_p p + e_n}.
$$

Substitute the detailed-balance relations $e_n/c_n = n_1$, $e_p/c_p = p_1$
(6.4):

$$
f_t = \frac{c_n n + c_p p_1}{c_n(n + n_1) + c_p(p + p_1)}.
$$

The net rate from (6.3):

$$
R = c_n n(1 - f_t) - c_n n_1 f_t = c_n n - c_n(n + n_1)f_t.
$$

Substitute $f_t$:

$$
R = c_n n - c_n(n + n_1)\cdot\frac{c_n n + c_p p_1}{c_n(n + n_1) + c_p(p + p_1)}.
$$

Common denominator:

$$
R = \frac{c_n n[c_n(n + n_1) + c_p(p + p_1)] - c_n(n + n_1)(c_n n + c_p p_1)}
        {c_n(n+n_1) + c_p(p+p_1)}.
$$

Expand the numerator:
$c_n^2 n(n+n_1) + c_nc_p n(p+p_1) - c_n^2 n(n+n_1) - c_nc_p p_1(n+n_1)$
$= c_nc_p[n(p+p_1) - p_1(n + n_1)]$
$= c_nc_p[np + np_1 - p_1 n - p_1 n_1]$
$= c_nc_p[np - n_1 p_1] = c_nc_p[np - n_i^2]$

(using $n_1 p_1 = n_i^2$). So

$$
R = \frac{c_nc_p[np - n_i^2]}{c_n(n+n_1) + c_p(p+p_1)}.
$$

Define $\tau_n = 1/(c_n N_t)$, $\tau_p = 1/(c_p N_t)$ (lifetimes).
Multiply numerator and denominator by $1/(c_nc_p N_t^2) = \tau_n\tau_p$:

$$
R = \frac{(np - n_i^2)/N_t^2}
        {(c_p^{-1}\tau_n^{-1})(n+n_1)/N_t + (c_n^{-1}\tau_p^{-1})(p+p_1)/N_t}.
$$

This simplifies via $1/c_n = \tau_n N_t$ etc. to give

$$
R = \frac{np - n_i^2}{\tau_p(n+n_1) + \tau_n(p+p_1)},
$$

matching (6.5). ✓

## B.3 — Depletion approximation with charge sheets

Goal: derive (7.8) for the depletion width $W$.

In the depletion region, the Poisson equation (in 1D) is

$$
\frac{d^2\psi}{dx^2} = -\frac{q}{\varepsilon}N_\mathrm{net}(x).
$$

For an abrupt step junction with $N_A$ on $-x_p < x < 0$ and $N_D$ on
$0 < x < x_n$, $N_\mathrm{net} = -N_A$ on the p-side and $+N_D$ on the
n-side. Integrate once:

$$
E^p(x) = -\frac{d\psi}{dx} = \frac{qN_A}{\varepsilon}(x + x_p)\;\;(-x_p < x < 0),
$$

$$
E^n(x) = \frac{qN_D}{\varepsilon}(x_n - x)\;\;(0 < x < x_n),
$$

with the boundary conditions $E(\pm x_n, \pm x_p) = 0$ (field vanishes
at the depletion edge). Continuity of $E$ at $x = 0$:

$$
\frac{qN_A x_p}{\varepsilon} = \frac{qN_D x_n}{\varepsilon}
\;\;\Rightarrow\;\; N_A x_p = N_D x_n.
\tag{B.2}
$$

This is charge balance.

Integrate $E^p$ from $-x_p$ to $0$:

$$
\psi(0) - \psi(-x_p) = -\int_{-x_p}^0 E^p\,dx
   = -\frac{qN_A}{\varepsilon}\int_{-x_p}^0 (x+x_p)\,dx
   = -\frac{qN_A x_p^2}{2\varepsilon}.
$$

Since $\psi(-x_p) = \psi_L^\mathrm{eq}$ (the bulk p-side equilibrium
value), $\psi(0) = \psi_L^\mathrm{eq} - qN_A x_p^2/(2\varepsilon)$.
But the convention is *positive* potential drop: take
$\psi_R^\mathrm{eq} - \psi_L^\mathrm{eq} = V_{bi} - V$ as the total
band bending.

$$
\psi_R^\mathrm{eq} - \psi_L^\mathrm{eq}
= \frac{qN_A x_p^2}{2\varepsilon} + \frac{qN_D x_n^2}{2\varepsilon}
= V_{bi} - V.
\tag{B.3}
$$

Combine (B.2) and (B.3) with $W = x_p + x_n$:

From (B.2): $x_p = (N_D/N_A)x_n$, so $x_n = (N_A/(N_A+N_D))W$ and
$x_p = (N_D/(N_A+N_D))W$.

Substitute into (B.3):
$N_A x_p^2 + N_D x_n^2 = N_A(N_D W/(N_A+N_D))^2 + N_D(N_A W/(N_A+N_D))^2$
$= W^2 N_A N_D(N_D + N_A)/(N_A+N_D)^2 = W^2 N_A N_D/(N_A+N_D)$.

So $V_{bi} - V = q W^2 N_AN_D/(2\varepsilon(N_A+N_D))$, giving

$$
W = \sqrt{\frac{2\varepsilon(V_{bi} - V)(N_A+N_D)}{qN_AN_D}}.
$$

(7.8) ✓

## B.4 — Slotboom transformation algebra

Goal: derive (11.6), $\mathbf{J}_n = -q\mu_n n\,\nabla\Phi_n$.

Start with (5.5):
$\mathbf{J}_n = q\mu_n n\,\mathbf{E} + qD_n\,\nabla n
= -q\mu_n n\,\nabla\psi + q\mu_n V_t\,\nabla n$ (using Einstein $D_n = \mu_n V_t$).

Slotboom: $n = n_i\exp((\psi - \Phi_n)/V_t)$, so

$$
\nabla n = n_i\exp((\psi - \Phi_n)/V_t)\cdot \frac{\nabla\psi - \nabla\Phi_n}{V_t}
        = \frac{n}{V_t}(\nabla\psi - \nabla\Phi_n).
$$

Substitute:

$$
\mathbf{J}_n = -q\mu_n n\,\nabla\psi + q\mu_n V_t\cdot\frac{n}{V_t}(\nabla\psi - \nabla\Phi_n)
            = -q\mu_n n\,\nabla\psi + q\mu_n n\,\nabla\psi - q\mu_n n\,\nabla\Phi_n
            = -q\mu_n n\,\nabla\Phi_n.
$$

The drift and diffusion terms cancelled. ✓

## B.5 — Dimensionless Poisson

Goal: derive (12.1).

Start with the dimensional Poisson (1.9):
$-\nabla\!\cdot\!(\varepsilon_0\varepsilon_r\,\nabla\psi) = q(p - n + N)$.

Substitute $\psi = V_t\hat\psi$, $n = C_0\hat n$, $p = C_0\hat p$,
$N = C_0\hat N$:

$$
-\nabla\!\cdot\!(\varepsilon_0\varepsilon_r V_t\,\nabla\hat\psi) = qC_0(\hat p - \hat n + \hat N).
$$

Divide both sides by $qC_0$:

$$
-\nabla\!\cdot\!\left(\frac{\varepsilon_0 V_t}{qC_0}\,\varepsilon_r\,\nabla\hat\psi\right)
= \hat p - \hat n + \hat N.
$$

Defining $L_D^2 = \varepsilon_0 V_t/(qC_0)$, the result is (12.1). ✓

The mesh stays in physical meters, so $\nabla$ is the physical gradient
(1/m). The coefficient $L_D^2 \cdot \varepsilon_r$ has units (m²)(dimensionless) =
m². Multiplying by $\nabla$ twice gives a dimensionless $\hat\psi$, as
required.

## B.6 — r-weighted axisymmetric Poisson

Goal: derive (15.4).

In cylindrical coordinates with $\theta$-independent fields, the
volume measure is $dV = r\,d\theta\,dr\,dz$. A 3D volume integral
$\int_\Omega f\,dV = 2\pi\int_{\Omega_\mathrm{merid}} f\,r\,dr\,dz$
where the integrand $f$ no longer depends on $\theta$.

The cylindrical divergence of a $\theta$-independent vector field
$\mathbf{F} = F_r\hat{\mathbf{e}}_r + F_z\hat{\mathbf{e}}_z$ is

$$
\nabla\!\cdot\!\mathbf{F} = \frac{1}{r}\frac{\partial(rF_r)}{\partial r} + \frac{\partial F_z}{\partial z}.
$$

The strong Poisson is $-\nabla\!\cdot\!(\varepsilon\nabla\psi) = \rho$.
Multiply by $v(r,z)$ and integrate against $r\,dr\,dz$:

$$
-\int\nabla\!\cdot\!(\varepsilon\nabla\psi)\cdot v\cdot r\,dr\,dz = \int\rho v r\,dr\,dz.
$$

Apply the divergence-theorem identity in cylindrical coordinates:
$\int_\Omega(\nabla\!\cdot\!\mathbf{F})v\,dV = -\int_\Omega \mathbf{F}\!\cdot\!\nabla v\,dV
+ \int_{\partial\Omega}\mathbf{F}\!\cdot\!\hat{\mathbf{n}}\,v\,dS$. With
$\mathbf{F} = \varepsilon\nabla\psi$:

$$
\int_\Omega\varepsilon\nabla\psi\!\cdot\!\nabla v\,dV - \int_{\partial\Omega}\varepsilon\nabla\psi\!\cdot\!\hat{\mathbf{n}}\,v\,dS
= \int_\Omega\rho v\,dV.
$$

Substituting $dV = r\,dr\,d\theta\,dz$ and dropping the $2\pi$:

$$
\int\varepsilon\nabla\psi\!\cdot\!\nabla v\,r\,dr\,dz - (\text{boundary terms})
= \int\rho v\,r\,dr\,dz.
$$

For homogeneous Neumann at $r = R$ and "natural" at $r = 0$, the
boundary terms drop. The result is (15.4). ✓

## B.7 — BDF2 truncation error

Goal: derive the leading truncation error of (17.5).

For $u(t)$ smooth, expand around $t^{n+1}$:

$$
u^n = u(t^{n+1} - dt) = u^{n+1} - dt\,u' + \frac{dt^2}{2}u'' - \frac{dt^3}{6}u''' + \frac{dt^4}{24}u'''' - ...
$$

$$
u^{n-1} = u^{n+1} - 2dt\,u' + 2dt^2\,u'' - \frac{4dt^3}{3}u''' + \frac{2dt^4}{3}u'''' - ...
$$

Plug into (17.5):
$\frac{1}{dt}\left[\frac{3}{2}u^{n+1} - 2u^n + \frac{1}{2}u^{n-1}\right]$
$= \frac{1}{dt}\bigg[\frac{3}{2}u^{n+1}$
$- 2(u^{n+1} - dt\,u' + \frac{dt^2}{2}u'' - \frac{dt^3}{6}u''' + ...)$
$+ \frac{1}{2}(u^{n+1} - 2dt\,u' + 2dt^2\,u'' - \frac{4dt^3}{3}u''' + ...)\bigg]$.

Collect coefficients:
$u^{n+1}$: $3/2 - 2 + 1/2 = 0$. ✓
$u'$: $0 - (-2dt) - dt = 2dt - dt = dt$. So $u'\cdot 1$ — first-order accuracy in $u'$ as expected.
$u''$: $0 - dt^2 + dt^2 = 0$. Second-order accuracy.
$u'''$: $0 - (-\frac{dt^3}{3}) - \frac{2dt^3}{3} = \frac{dt^3}{3} - \frac{2dt^3}{3} = -\frac{dt^3}{3}$.
Divided by $dt$: $-dt^2 u'''/3$. So the truncation is $-(dt^2/3) u'''$
to leading order, *second-order accuracy*. ✓ (BDF2 is $O(dt^2)$.)

(Different references quote $-(2/9)dt^2 u'''$ depending on
normalization; the key point is the $dt^2$ scaling.)

## B.8 — AC linearization with phasor algebra

Goal: derive (18.2) and the real 2×2 block (18.4).

Let $u(t) = u_0 + \delta u\,e^{j\omega t}$, $V(t) = V_\mathrm{DC} + \delta V\,e^{j\omega t}$.
The transient residual is

$$
F_\mathrm{trans}(u; V) = F_\mathrm{ss}(u; V) + M\,\dot u = 0.
$$

Linearize around $(u_0, V_\mathrm{DC})$:

$$
F_\mathrm{ss}(u_0 + \delta u; V_\mathrm{DC} + \delta V)
\approx F_\mathrm{ss}(u_0; V_\mathrm{DC}) + J\delta u + \frac{\partial F}{\partial V}\delta V
= 0 + J\delta u + \frac{\partial F}{\partial V}\delta V.
$$

The time derivative: $\dot u = j\omega\delta u\,e^{j\omega t}$. So

$$
J\delta u + \frac{\partial F}{\partial V}\delta V + j\omega M\delta u = 0,
$$

giving (18.2):

$$
(J + j\omega M)\delta u = -\frac{\partial F}{\partial V}\delta V.
$$

Real-block reformulation: $\delta u = x + jy$, RHS $= b_R + jb_I$:
$(J + j\omega M)(x + jy) = b_R + jb_I$
$Jx + j\omega Mx + jJy - \omega My = b_R + jb_I$
$(Jx - \omega My) + j(\omega Mx + Jy) = b_R + jb_I$.

Equate real and imaginary:
$Jx - \omega My = b_R$
$\omega Mx + Jy = b_I$.

In matrix form:
$\begin{bmatrix} J & -\omega M \\ \omega M & J\end{bmatrix}\begin{bmatrix}x\\y\end{bmatrix}
= \begin{bmatrix}b_R\\b_I\end{bmatrix}$.

Equation (18.4). ✓

## See also

- Each derivation has a chapter cross-reference; see the chapter for the
  motivation and consequence of the result.
- Appendix C — symbol/glossary for any unfamiliar variable.
- For implementation details, see the Code Map sections in each chapter.

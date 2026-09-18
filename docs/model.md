# Model reference

## Transport

For mass fraction `m`, the maintained solvers use

- `rho = rho0 * (1 + rho_mACoeff*m)`;
- `mu = mu0 * (1 + mu_mACoeff*m)`;
- `D = max(D_AB_Coeff*(1 - D_AB_mACoeff*m), D_AB_Min)`.

The 2012 density expression `997.1 + 694*m` therefore requires
`rho_mACoeff = 694/997.1`, not 694. Osmotic pressure is `80510000*m` Pa.
The printed 2012 viscosity exponent differs from the physical SI value and from
the later 2016 paper: use `mu0 = 0.00089` Pa s, as in the original repository.
The empirical NaCl properties apply near 25 C and mass fractions up to 0.09.

## Membrane equations

Let `j >= 0` be feed-side water velocity normal to the membrane, `f` and `d` the
feed/draw surface mass fractions, and `P` the osmotic-pressure coefficient. The
finite-permeability model (2012/2016 Eq. 4) is evaluated using the equivalent
residual

```
F(j) = j + B + A*P*f - (B + A*P*d)*exp(-j*K).
```

This avoids the logarithmic singularity at zero concentration and the division
by zero at `K=0`. `F` is monotone for the supported physical parameters, with a
root in `[0, A*P*(d-f)]`. Safeguarded Newton iterations retain this bracket and
fall back to bisection when necessary. The absolute velocity tolerance is
1e-16 m/s plus a relative tolerance of 1e-12. Since `F′ >= 1`, the residual also
bounds the flux error. Roots are reused only for exactly unchanged surface
concentrations; velocity, density scaling, and slip updates still run. For
`eq simple`, omit `B` from this residual. The two models coincide at `B=0`; at `K=0`, both give the no-ICP limit.

Water mass is conserved by scaling the draw-side velocity by `rho_feed/rho_draw`.
The salt-flux relation uses the papers' dilute-solution approximation
`phi_osmotic = P/1000` Pa m3/kg. The paired reverse salt flux has the same magnitude
on both sides, based on the feed-side water velocity. For outward normal `n`,
transport leaving a fluid compartment is

```
q_s = rho*m*(U dot n) - rho*D*(dm/dn).
```

Salt leaves the draw compartment and enters the feed compartment.

## Boundary and solver behavior

Membrane faces are paired by position, independently of patch-face ordering.
Each pair must have equal area and opposite normals and remain on one MPI rank.
Normal velocities and slip use dot-product projections.

SIMPLE uses a momentum predictor and consistent pressure-gradient discretization
in the predictor and corrector. Fixed membrane and inlet velocities are preserved
in the pressure predictor; the velocity corrector uses relaxed pressure.
Steady and transient solvers use the same concentration-dependent viscosity law.

`potentialSalt` evaluates flow-rate inlet conditions before initialization.
Use `-writep` to write physical pressure.

A valid scalar boundary layer still requires sufficient resolution. The explicit
RO Robin condition can become singular when `U_n*R*delta/D` approaches one;
refine the wall-normal grid instead of interpreting an unstable run as physical.
FO rejects materially negative surface mass fractions; values within 1e-12 below
zero are treated as roundoff only when evaluating its flux.

## Steady convergence and conservation

`membraneSimpleControl` extends SIMPLE's residual checks with the actual discrete
salt flux, including non-orthogonal diffusion corrections. The absolute net
boundary salt flux is divided by one-sided membrane salt transfer. The default
limit is 0.01 (1%), configured as `maxRelativeSaltImbalance` in `SIMPLE`.
This catches premature convergence when a small error relative to bulk inlet
salt is still large relative to the much smaller membrane transfer. The integral
check supplements the residual limits; it does not replace them.

The transient solver reports boundary fluxes too, but transient salt storage
must be included before interpreting those fluxes as a conservation error.

# Consistent pressure initialization

`simpleSaltTransport` provides an optional consistent pressure correction for
steady initialization. It is disabled by default and supports CPU and GPU
execution. The [GPU execution guide](gpu.md) gives the chamber startup and
original-equation convergence commands.

## Pressure correction

Let `rAU` be the inverse relaxed momentum diagonal and define
`rAtU = 1 / (1/rAU - H1)`, using OpenFOAM's `UEqn().H1()`. With the current
pressure and density held fixed, the predictor receives
`(rAtU - rAU) * (grad(p) + gravityDensityForce)`. The velocity corrector and
pressure Laplacian use `rAtU`. The face predictor receives the flux of the
discrete pressure matrix with coefficient `rho * (rAtU - rAU)`, including
boundary and non-orthogonal terms. This preserves the original same-pressure
cell velocity and pressure-corrected face flux under the conditions below.

## Supported discretization

Use a static mesh, linear pressure boundary conditions and compatible pressure
Laplacians with linear coefficient interpolation. The supported chamber scheme
is `Gauss linear corrected`. Coefficient-specific Laplacian overrides must use
compatible discretizations. Nonlinear coefficient interpolation, including
harmonic interpolation, is outside the validated scope.

The implementation accepts `fixedValue`, `zeroGradient`, `empty`, `processor`,
`symmetryPlane` and `symmetry` pressure patches and rejects other patch types.
It rejects non-positive or non-finite modified momentum diagonals. These guards
do not verify every discretization assumption.

## Identity verification

Enable the optional identity check in `system/fvSolution`:

```foam
SIMPLE
{
    // Retain the case's other SIMPLE settings.
    consistent      true;
    consistentAudit true;
}
```

The audit compares original and modified same-pressure cell velocities and
checks `newPressureFlux - oldPressureFlux - deltaPressureFlux` on internal and
boundary faces. Each relative defect must be at most `1e-10`; a failed identity
stops the run. Disable `consistentAudit` after qualifying the discretization,
or include its cost consistently in CPU/GPU timings.

## Final convergence

Consistent startup can use upwind salt convection and different relaxation
factors. Restore the intended `fvSchemes`, relaxation factors, linear targets
and iteration caps at the handoff, retain the written fields and set
`consistent false`. Final convergence applies to those restored equations.

`SIMPLE.validationConsecutiveChecks` requires consecutive successful original
outer residual and salt-balance checks. Its default is `1`; valid values are
1 to 100,000. A failed residual or salt check resets the counter. Physical mass
balance, water-flux stability, concentration bounds and local-field accuracy
require independent checks, including when the startup identity passes.

The initialization method does not establish a unique steady solution.
Validate final local fields and integral fluxes for the intended problem.

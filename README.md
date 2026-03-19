# Brinkman Penalization — solid disk entrained by channel flow

Basilisk simulations of a non-rotating solid disk freely advected by a
2-D plug-flow channel, using the **Brinkman Penalization Method (BPM)**.

---

## Physical problem

A rectangular channel (height H = 2, length H = 2) carries a uniform
inlet flow at mean velocity `Umean`.  A solid disk of radius `Rc` is
placed near the inlet and released from rest.  The flow accelerates the
disk; the disk does not rotate.  Both fluid and solid have density
ρ = 1.

The Reynolds number is defined as Re = Umean · H / ν.

---

## Method — Brinkman Penalization (BPM)

The solid's no-slip condition is enforced *without* an explicit
body-fitted mesh.  Instead, a penalization body force is added to the
Navier–Stokes momentum equation everywhere inside the solid:

```
f_B = −χ/η · (u − u_s)
```

| Symbol | Meaning |
|--------|---------|
| χ      | Solid volume fraction (0 outside, 1 inside, computed via `fraction()`) |
| η      | Penalization parameter ≈ h²/ν (vanishing permeability at grid scale) |
| u_s    | Prescribed solid velocity (Ucx, Ucy) |

The reaction force on the disk (Newton's 3rd law) is integrated to
advance the disk centroid with a forward-Euler scheme.

---

## Scripts

### `binkmanSimple.c` — **recommended, best agreement with COMSOL**

The simplest and most accurate implementation.

- Uses the standard Basilisk solver `navier-stokes/centered.h` unchanged.
- The BPM force is applied as an **explicit face-centered acceleration**
  in the `acceleration` event — no modifications to the viscosity solver.
- Fluid and solid share the same density and kinematic viscosity, so no
  separate material-property fields are needed.
- Uses wavelet-based adaptive mesh refinement (AMR) on χ and u.

**Key parameters** (top of file):

| Variable | Default | Description |
|----------|---------|-------------|
| `Umean`  | 0.1 | Inlet plug-flow velocity |
| `Reynolds` | 100 | Re = Umean·H/ν |
| `TEND`   | 5.0 | End time |
| `DTOUT`  | 0.5 | Output snapshot interval |
| `LEVELINI` | 7 | Initial grid level (128² cells) |
| `maxlevel` | 9 | Max AMR level (512² equivalent) |
| `Rc`     | 0.125 | Disk radius |
| `xc, yc` | 0.25, 1.0 | Initial disk centre |

---

### `binkmanViscosity.c` — alternative, BPM inside viscosity solver

A more sophisticated variant where the penalization term is incorporated
directly into the **implicit multigrid viscosity solve**.

- Uses custom headers `centeredBinkman.h` and `viscosityBinkman.h`
  (modified copies of Basilisk's standard solvers).
- Variable-density and variable-viscosity Eulerian fields are computed
  as linear mixtures of fluid and solid values weighted by χ.
- Inlet velocity is ramped from 0 to `Umean` over the first 0.1 s.
- AMR is disabled by default.

> **Note:** `binkmanSimple.c` matches COMSOL reference results better
> than this variant.  `binkmanViscosity.c` is kept for methodological
> comparison.

---

### `centeredBinkman.h` — modified centered N-S solver

Drop-in replacement for `navier-stokes/centered.h`.  The only
difference is that it includes `viscosityBinkman.h` instead of
`viscosity.h`, routing the BPM fields into the implicit viscosity step.

Used only by `binkmanViscosity.c`.

---

### `viscosityBinkman.h` — modified implicit viscosity solver

Extends Basilisk's `viscosity.h` to include the BPM term

```
Kb = χ · dt / η
```

on the diagonal of the implicit viscous system.  Because Kb ≫ 1
inside the solid, the multigrid solve drives u → u_s there.

Declares three global fields that the user script must update each step:
- `bpm_chi[]` — solid volume fraction
- `bpm_eta[]` — penalization parameter
- `bpm_us[]`  — prescribed solid velocity

Used only by `binkmanViscosity.c` (via `centeredBinkman.h`).

---

## Building

The project uses the standard Basilisk build system.  Make sure the
`BASILISK` environment variable points to the Basilisk source tree.

```bash
# compile and run binkmanSimple
make binkmanSimple.tst

# or compile manually
qcc -O2 -disable-dimensions -DJACOBI=1 binkmanSimple.c -o binkmanSimple -lm
./binkmanSimple > log
```

The `Makefile` sets `-O2 -disable-dimensions -DJACOBI=1` by default
(Jacobi relaxation is preferred for the multigrid solver).

---

## Outputs

Both scripts write to the current working directory:

| File | Content |
|------|---------|
| `log` (stdout) | Time series: step, t, dt, xc, yc, Ucx, Ucy, Fhx, Fhy, mgp.i, mgu.i |
| `ux.png` | x-velocity field with disk outline (overwritten every `DTOUT`) |
| `omega.png` | Vorticity field with streamlines and disk outline |

Log column reference:

```
# 1:step 2:t 3:dt 4:xc 5:yc 6:Ucx 7:Ucy 8:Fhx 9:Fhy 10:mgp.i 11:mgu.i
```

---

## Dependencies

- [Basilisk](http://basilisk.fr) C flow solver
- Standard Basilisk headers: `navier-stokes/centered.h`, `fractions.h`, `view.h`
- Custom headers (this directory): `centeredBinkman.h`, `viscosityBinkman.h`

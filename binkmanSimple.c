/**
# Plug-flow channel with a freely advected solid disk — simple Brinkman penalization

## Physical setup

A 2-D channel of height H = 2 and length H = 2 carries a uniform
(plug) flow at mean velocity Umean.  A solid disk of radius Rc,
centred initially at (xc, yc), is released from rest and entrained by
the flow.  The disk **does not rotate**.

Reynolds number is defined as Re = Umean * H / nu_f.

## Numerical method — Brinkman Penalization Method (BPM)

The no-slip condition inside the solid is enforced weakly by adding a
volume penalization body force to the Navier-Stokes momentum equation:

    f_B = -chi/eta * (u - us)

where:
- chi   : volume fraction of the solid (0 outside, 0–1 in interfacial
          cells, 1 fully inside), computed via fraction().
- eta   : penalization parameter, eta = (h_min)^2 / nu_f, chosen so
          that the penalized region behaves as a porous medium with
          vanishing permeability at the grid scale.
- us    : prescribed solid velocity (Ucx, Ucy) — zero rotation assumed.

The force is added as an **explicit face-centered acceleration** in the
`acceleration` event, using face-averaged chi and eta values.  This
avoids modifying the viscosity solver and keeps the implementation
minimal.

## Material properties

Fluid and solid share the same density (rho = 1) and kinematic
viscosity (nu_f), so no separate material-property fields are needed.
The standard `navier-stokes/centered.h` solver is used unchanged.

## Solid body dynamics

At each timestep the reaction force on the disk (Newton's 3rd law
applied to f_B) is integrated with a forward-Euler scheme to advance
the disk centroid velocity (Ucx, Ucy) and position (xc, yc).

    Fhx = integral over solid of chi/eta * (u.x - Ucx) dV
    Ucx += dt * Fhx / mp
    xc  += dt * Ucx

where mp = pi * Rc^2 is the disk mass per unit depth (rho_s = 1).

## Mesh adaptivity

The quad-tree mesh is refined around the disk boundary (bpm_chi) and
velocity gradients (u) using wavelet-based adaptation between minlevel
and maxlevel.

## Outputs

- **stdout / log**: one line per timestep with columns
  step, t, dt, xc, yc, Ucx, Ucy, Fhx, Fhy, mgp.i, mgu.i
- **ux.png**   : x-velocity field with disk outline (updated every DTOUT)
- **omega.png**: vorticity field with streamlines and disk outline

## Validation

This simplified formulation (BPM force as explicit acceleration,
uniform material properties) agrees well with a reference COMSOL model
and outperforms the more complex `binkmanViscosity.c` variant, which
incorporates the penalization term directly inside the implicit
viscosity multigrid solver.

## Parameters (edit at the top of the file)

| Variable   | Default | Meaning                                   |
|------------|---------|-------------------------------------------|
| Umean      | 0.1     | Inlet plug-flow velocity                  |
| Reynolds   | 100.0   | Reynolds number Re = Umean*H/nu           |
| TEND       | 5.0     | End time                                  |
| DTOUT      | 0.5     | Output interval                           |
| LEVELINI   | 7       | Initial uniform grid level (128^2 cells)  |
| maxlevel   | 9       | Maximum AMR level (512^2 equivalent)      |
| minlevel   | 3       | Minimum AMR level                         |
| Rc         | 0.125   | Disk radius                               |
| xc, yc     | 0.25, 1 | Initial disk centre                       |
*/

//#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "fractions.h"
#include "view.h"

#define H  2.0
#define PI 3.14159265358979323846

/* Flow parameters */
double Umean   = 0.1;
double Reynolds = 100.0;
double TEND    = 5.0;
double DTOUT    = 0.5;
int    LEVELINI = 7;
int    maxlevel = 9;
int    minlevel = 3;

/* Disk (rho_s = rho_f = 1) */
double Rc  = 0.125;
double xc  = 0.25, yc = 1.0;
double Ucx = 0.0,  Ucy = 0.0;
double Fhx = 0.0,  Fhy = 0.0;
double mp  = 0.0;          /* disk mass per unit depth */

double nu_f = 0.0;

face vector muv[];

/* ------------------------------------------------------------------ */
/* Boundary conditions                                                  */
/* ------------------------------------------------------------------ */

u.n[left]   = dirichlet(Umean);
u.t[left]   = dirichlet(0.);
p[left]     = neumann(0.);
pf[left]    = neumann(0.);

u.n[right]  = neumann(0.);
u.t[right]  = neumann(0.);
p[right]    = dirichlet(0.);
pf[right]   = dirichlet(0.);

u.n[top]    = dirichlet(0.);
u.t[top]    = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);

psi[bottom] = dirichlet(0.);
psi[top]    = dirichlet(-2.*Umean);
psi[left]   = dirichlet(-Umean*y);
psi[right]  = dirichlet(-Umean*y);


/* ------------------------------------------------------------------ */
/* BPM helpers                                                          */
/* ------------------------------------------------------------------ */

/* Brinkman fields */
scalar bpm_chi[];
scalar bpm_eta[];
vector bpm_us[];
scalar psi[];

static inline void update_bpm_fields()
{
  fraction(bpm_chi, sq(Rc) - sq(x - xc) - sq(y - yc));

  /* eta = (m h_min)^2 / nu, uniform inside solid */
  double eta = sq(2. * L0 / (1 << maxlevel)) / (nu_f + SEPS);
  foreach() {
    bpm_us.x[] = Ucx;
    bpm_us.y[] = Ucy;
    bpm_eta[]  = (bpm_chi[] > 0.) ? eta : 0.;
  }
}

/* ------------------------------------------------------------------ */
/* Main                                                                 */
/* ------------------------------------------------------------------ */

int main()
{
  L0 = H;
  origin(0., 0.);
  init_grid(1 << LEVELINI);

  nu_f = Umean * H / Reynolds;
  mp   = PI * sq(Rc);       /* rho_s = 1 */

  DT = 1.e-2;
  mu = muv;
  run();
}

/* ------------------------------------------------------------------ */
/* Events                                                               */
/* ------------------------------------------------------------------ */

event properties(i++)
{
  foreach_face()
    muv.x[] = fm.x[] * nu_f;
  update_bpm_fields();
}

event defaults (i = 0) {
	if (is_constant (a.x))
		a = new face vector;

}

event init(t = 0)
{
  update_bpm_fields();
  foreach() {
    u.x[] = 0.;
    u.y[] = 0.;
  }
}

/**
Compute the BPM reaction force on the disk by integrating
  F_h = integral( chi/eta * (u - Us) dV )
over all cells with chi > 0, then update disk kinematics with a
forward-Euler step (Newton's 3rd law: the fluid exerts -F_h on the
solid, so the solid acceleration is F_h/mp). */

event end_timestep(i++)
{
  Fhx = Fhy = 0.;
  foreach(reduction(+:Fhx) reduction(+:Fhy)) {
    if (bpm_chi[] > 1e-12) {
      double coeff = bpm_chi[] * dv() / (bpm_eta[] + SEPS);
      Fhx += coeff * (u.x[] - Ucx);
      Fhy += coeff * (u.y[] - Ucy);
    }
  }

  /* Forward-Euler rigid body update */
  xc  += dt * Ucx;
  yc  += dt * Ucy;
  Ucx += dt * Fhx / mp;
  Ucy += dt * Fhy / mp;
}

event logfile(i++; t <= TEND)
{
  if (i == 0)
    fprintf(stdout,
            "# 1:step 2:t 3:dt 4:xc 5:yc 6:Ucx 7:Ucy"
            " 8:Fhx 9:Fhy 10:mgp.i 11:mgu.i\n");
  fprintf(stdout, "%d %g %g %g %g %g %g %g %g %d %d\n",
          i, t, dt, xc, yc, Ucx, Ucy, Fhx, Fhy, mgp.i, mgu.i);
}

event snapshot(t = 0.; t += DTOUT; t <= TEND)
{
  char label[80];

  clear();
  view(tx = -0.5, ty = -0.5, bg = {1, 1, 1});
  box();
	sprintf(label, "Ux t:%3.2f", t);
  cells();
	squares("u.x", min = 0, max = 1.5, cbar = true, pos = {0.5, -0.2}, levels = 10, mid = true, format = "%6.4f", label = label);
  //squares("u.x", linear = true, min = 0., max = Umean);
  isoline("bpm_chi", 0.5, lw = 2);
  //vectors("u", scale = 0.5);
  save("ux.png");

  clear();
  scalar omega[];
  vorticity (u, omega);
  poisson ( psi, omega);
	sprintf(label, "Ux t:%3.2f", t);
  squares ("omega", linear = true, cbar = true, pos = {0.5, -0.2}, levels = 10, mid = true, format = "%6.4f", label = label);
  isoline ("psi", n = 50, lw = 1);
  isoline ("bpm_chi", 0.5 , lw = 3);
  cells();
  //draw_vof ("bpm_chi", lc = {0,0,0}, lw = 2);
  save ("omega.png");
}

event acceleration (i++)
{
  face vector av = a;
  foreach_face(x) {
    double chif = (bpm_chi[] + bpm_chi[-1])/2.;
    if (chif > 1e-6)
      av.x[] -= chif / ((bpm_eta[] + bpm_eta[-1])/2. + SEPS) * (u.x[] - (bpm_us.x[] + bpm_us.x[-1])/2.);
  }
  foreach_face(y) {
    double chif = (bpm_chi[0,0] + bpm_chi[0,-1])/2.;
    if (chif > 1e-6)
      av.y[] -= chif / ((bpm_eta[0,0] + bpm_eta[0,-1])/2. + SEPS) * (u.y[] - (bpm_us.y[0,0] + bpm_us.y[0,-1])/2.);
  }
}



#if 1
event adapt(i++) {
  adapt_wavelet((scalar *){bpm_chi,u},
                (double []){1e-3, 1e-3, 1e-3},
                maxlevel, minlevel);
}
#endif
/**
# Plug-flow channel with a freely advected solid disk — Brinkman penalization via implicit viscosity

## Physical setup

Same geometry as `binkmanSimple.c`: a 2-D channel (H × H) with a
uniform inlet flow at Umean entraining a non-rotating solid disk of
radius Rc.

## Numerical method — BPM coupled to the implicit viscosity solver

Unlike `binkmanSimple.c`, this variant incorporates the Brinkman
penalization term directly into the **implicit multigrid viscosity
solve** (see `viscosityBinkman.h`).  The penalization coefficient

    Kb = chi * dt / eta,    eta = (m_brink * Delta)^2 / nu_f

is added to the diagonal of the viscous system, driving the velocity
inside the solid toward the prescribed solid velocity us = (Ucx, Ucy)
in the same implicit sweep that handles viscous diffusion.  This
coupling is tighter than the explicit-acceleration approach used in
`binkmanSimple.c`.

The custom solver chain is:
    centeredBinkman.h  →  (replaces navier-stokes/centered.h)
    viscosityBinkman.h →  (replaces viscosity.h; contains BPM terms)

## Material properties

Variable-density / variable-viscosity fields (rhov, muc, alphav, muv)
are computed each timestep as linear mixtures of fluid and solid values
weighted by bpm_chi.  Fluid and solid are given the same kinematic
viscosity (mu_s = rho_s * nu_f) but separate dynamic viscosities are
maintained for generality.

## Solid body dynamics

Identical forward-Euler Newton integration as in `binkmanSimple.c`,
but the penalization force density is computed as

    f = rho * (u - us) / eta

(includes the local density, unlike the simpler variant).

## Mesh adaptivity

Adaptive refinement is **disabled by default** (`#if 0` block).

## Outputs

- **stdout / log**: step, t, dt, xc, yc, Ucx, Ucy, Fhx, Fhy, mgp.i, mgu.i
- **ux.png**   : x-velocity field with disk boundary (updated every DTOUT)
- **omega.png**: vorticity field with streamlines and disk boundary

## Comparison with binkmanSimple.c

`binkmanSimple.c` (BPM as explicit acceleration, uniform properties)
matches COMSOL reference results better than this more complex variant.
This file is kept for methodological comparison.

## Parameters (edit at the top of the file)

| Variable   | Default | Meaning                                      |
|------------|---------|----------------------------------------------|
| Umean      | 0.1     | Inlet velocity (ramp over first 0.1 s)       |
| Reynolds   | 100.0   | Re = Umean*H/nu                              |
| TEND       | 10.0    | End time                                     |
| DTOUT      | 0.1     | Output interval                              |
| LEVELINI   | 7       | Initial grid level                           |
| maxlevel   | 9       | Maximum AMR level                            |
| m_brink    | 1.0     | BPM coefficient: eta = (m_brink*Delta)^2/nu  |
| Rc         | 0.125   | Disk radius                                  |
| xc, yc     | 0.25, 1 | Initial disk centre                          |
*/

#include "grid/multigrid.h"
#include "centeredBinkman.h"
#include "fractions.h"
#include "view.h"

#define H 2.0
#define PI 3.14159265358979323846

/* Channel flow */
double Umean = 0.1;
double Reynolds = 100.0;
double TEND = 10.0;
double DTOUT = 0.1;
int LEVELINI = 7;
int maxlevel = 9;
int minlevel = 3;

double rho_f = 1.0;
double rho_s = 1.0;
double m_brink = 1.; /* eta = (m h)^2/nu */

/* Disk */
double Rc = 0.125;
double xc = 0.25, yc = 1.0;
double Ucx = 0.0, Ucy = 0.0;
double Fhx = 0.0, Fhy = 0.0;
double mp = 0.0;

/* Eulerian material fields */
scalar rhov[];
scalar muc[];
face vector alphav[];
face vector muv[];
scalar psi[];

double nu_f = 0.0, mu_f = 0.0, mu_s = 0.0;

#define UINLET(y) (6.0*Umean*((y)/H)*(1.0 - (y)/H))

static inline double uInlet ( double y, double t, double Um, double h)
{
  return Um*(t<0.1 ? t/0.1 : 1.0);
  //return 6.*Um*(y/h)*(1-y/h)*(t<1.0 ? t : 1.0);
}

u.n[left]  = dirichlet(uInlet(y, t, Umean, H));
u.t[left]  = dirichlet(0.0);
p[left]    = neumann(0.0);
pf[left]   = neumann(0.0);

u.n[right] = neumann(0.0);
u.t[right] = neumann(0.0);
p[right]   = dirichlet(0.0);
pf[right]  = dirichlet(0.0);

u.n[top] = dirichlet(0.0);
u.t[top] = neumann(0.0);

u.n[bottom] = dirichlet(0.0);
u.t[bottom] = neumann(0.0);

psi[bottom] = dirichlet(0.);
psi[top]    = dirichlet(-2.*Umean);
psi[left]   = dirichlet(-Umean*y);
psi[right]  = dirichlet(-Umean*y);

static inline void update_disk_mask()
{
  fraction (bpm_chi, sq(Rc) - sq(x - xc) - sq(y - yc));
}

static inline void update_material_fields()
{
  update_disk_mask();

  foreach() {
    bpm_us.x[] = Ucx;
    bpm_us.y[] = Ucy;
    bpm_eta[]  = (bpm_chi[] > 0.) ? sq(m_brink*Delta)/(nu_f + SEPS) : 0.;

    rhov[] = (1.0 - bpm_chi[])*rho_f + bpm_chi[]*rho_s;
    muc[]  = (1.0 - bpm_chi[])*mu_f  + bpm_chi[]*mu_s;
  }
  //boundary ((scalar *) {bpm_chi, bpm_eta, bpm_us, rhov, muc});

  foreach_face() {
    double rf = (rhov[] + rhov[-1])/2.;
    double muf = (muc[] + muc[-1])/2.;
    alphav.x[] = fm.x[]/(rf + SEPS);
    muv.x[]    = fm.x[]*muf;
  }
}

int main()
{
  L0 = H;
  origin (0.0, 0.0);
  init_grid (1 << LEVELINI);

  nu_f = Umean*H/Reynolds;
  mu_f = rho_f*nu_f;
  mu_s = rho_s*nu_f; /* same kinematic viscosity as the liquid */
  //mu_s = 2.*mu_f;
  mp   = rho_s*PI*sq(Rc);

  DT = 1.e-4;
  mu = muv;
  alpha = alphav;
  rho = rhov;
  run();
}

event init (t = 0)
{
  Ucx = 0.0;
  Ucy = 0.0;

  update_material_fields();

  foreach() {
    u.x[] = 0.0;
    u.y[] = 0.0;
  }
}

event logfile (i++; t <= TEND)
{
  if (i == 0)
    fprintf (stdout,
             "# 1:step 2:t 3:dt 4:xc 5:yc 6:Ucx 7:Ucy 8:Fhx 9:Fhy 10:mgp.i 11:mgu.i\n");
  fprintf (stdout,
           "%d %g %g %g %g %g %g %g %g %d %d\n",
           i, t, dt, xc, yc, Ucx, Ucy, Fhx, Fhy, mgp.i, mgu.i);
}

event properties (i++)
{
  update_material_fields();
}

event dynamics (i++)
{
  double Fx = 0., Fy = 0.;
  foreach(reduction(+:Fx) reduction(+:Fy)) {
    if (bpm_eta[] > 0.) {
       // Force density on solid = rho * (u - u_s) / eta
       double f_dens_x = rhov[] * (u.x[] - bpm_us.x[]) / bpm_eta[];
       double f_dens_y = rhov[] * (u.y[] - bpm_us.y[]) / bpm_eta[];
       Fx += f_dens_x * dv();
       Fy += f_dens_y * dv();
    }
  }
  Fhx = Fx;
  Fhy = Fy;

  double dt_p = dt;
  Ucx += (Fhx / mp) * dt_p;
  Ucy += (Fhy / mp) * dt_p;

  xc += Ucx * dt_p;
  yc += Ucy * dt_p;
}

event snapshot (t=0; t += DTOUT; t <= TEND)
{
  char name[80];
  char label[80];
  sprintf (name, "snapshot-%0.4f.png", t);

  view (tx = -0.5, ty = -0.5,
        samples = 4,
        bg = {1,1,1});
  box();


  clear();
  box();
  sprintf (name, "ux-%0.4f.png", t);
	sprintf(label, "Ux t:%3.2f", t);
	squares("u.x", min = 0, max = 0.15, cbar = true, pos = {0.5, -0.2}, levels = 10, mid = true, format = "%6.4f", label = label);
  draw_vof ("bpm_chi", lc = {0,0,0}, lw = 2);
  //save (name);
  save ("ux.png");

  clear();
  scalar omega[];
  vorticity (u, omega);
  poisson ( psi, omega);
	sprintf(label, "Ux t:%3.2f", t);
  squares ("omega", linear = true, min=statsf(omega).min, max=statsf(omega).max, cbar = true, pos = {0.5, -0.2}, levels = 10, mid = true, format = "%6.4f", label = label);
  isoline ("psi", n = 50, lw = 1);
  isoline ("bpm_chi", 0.5 , lw = 3);
  //cells();
  //draw_vof ("bpm_chi", lc = {0,0,0}, lw = 2);
  save ("omega.png");
}

#if 0
event adapt (i++)
{
  adapt_wavelet ((scalar *){u},
                 (double []){1e-3, 1e-3},
                 maxlevel, minlevel);
}
#endif

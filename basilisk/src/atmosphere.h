#include "utils.h"

face vector u[], un[];
scalar h[], hn[], zb[];

// Default parameters
// Coriolis parameter
double F0 = 1.;
// acceleration of gravity
double G = 1.;
// Viscosity
double NU = 0.;

/**
 * @brief Computes the centered finite difference advection of a scalar field.
 *
 * Calculates the advection term for the scalar field @p f using the velocity field @p u,
 * applying a centered (symmetric) finite difference scheme. The result is stored in @p df.
 */
trace
void advection_centered (scalar f, vector u, scalar df)
{
  foreach()
    df[] = ((f[] + f[-1,0])*u.x[] - 
	    (f[] + f[1,0])*u.x[1,0] +
	    (f[] + f[0,-1])*u.y[] - 
	    (f[] + f[0,1])*u.y[0,1])/(2.*Delta);
}

/**
 * @brief Computes the upwind finite difference advection of a scalar field.
 *
 * Calculates the advection term for the scalar field `f` using the upwind scheme based on the velocity field `u`, and stores the result in `df`. The upwind direction is chosen according to the sign of each velocity component to enhance numerical stability.
 */
trace
void advection_upwind (scalar f, vector u, scalar df)
{
  foreach()
    df[] = ((u.x[] < 0. ? f[] : f[-1,0])*u.x[] - 
	    (u.x[1,0] > 0. ? f[] : f[1,0])*u.x[1,0] +
	    (u.y[] < 0. ? f[] : f[0,-1])*u.y[] - 
	    (u.y[0,1] > 0. ? f[] : f[0,1])*u.y[0,1])/Delta;
}

trace
double /**
 * @brief Computes the maximum stable timestep for the shallow water system based on CFL conditions.
 *
 * Evaluates the timestep constraints imposed by gravity wave speed and local velocity magnitudes across all grid cells, returning the largest permissible timestep that maintains numerical stability.
 *
 * @return double The computed stable timestep.
 */
timestep (void)
{
  double dtmax = DT/CFL;
  dtmax *= dtmax;
  foreach(reduction(min:dtmax)) {
    Delta *= Delta;
    if (h[] > 0.) {
      double dt = Delta/(G*h[]);
      if (dt < dtmax) dtmax = dt;
    }
    foreach_dimension()
      if (u.x[] != 0.) {
	double dt = Delta/sq(u.x[]);
	if (dt < dtmax) dtmax = dt;
      }
  }
  return sqrt (dtmax)*CFL;
}

trace
void momentum (vector u, scalar h, vector du)
{
  scalar ke[];
  vertex scalar psi[];
  scalar dux[], dvy[];
  vector d;
  d.x = dux; d.y = dvy;

  foreach() {
#if 1
    ke[] = (sq(u.x[] + u.x[1,0]) + sq(u.y[] + u.y[0,1]))/8.;
#else
    double uc = u.x[]*u.x[] + u.x[1,0]*u.x[1,0];
    double vc = u.y[]*u.y[] + u.y[0,1]*u.y[0,1];
    ke[] = (uc + vc)/4.;
#endif
    foreach_dimension()
      d.x[] = (u.x[1,0] - u.x[])/Delta;
  }
  foreach_vertex()
    psi[] = (u.y[] - u.y[-1,0] + u.x[0,-1] - u.x[])/Delta;  

  coord f = {1.,-1.};
  foreach_face()
    du.x[] = 
      - (G*(h[] + zb[]) + ke[] - G*(h[-1,0] + zb[-1,0]) - ke[-1,0])/Delta
      + f.x*(((psi[] + psi[0,1])/2. + F0)*
	     (u.y[] + u.y[0,1] + u.y[-1,0] + u.y[-1,1])/4.)
      + NU*(u.x[0,1] + u.x[0,-1] - 2.*u.x[])/sq(Delta)
      + NU*(d.x[] - d.x[-1,0])/Delta;
}

/**
 * @brief Computes time derivatives of velocity and height fields for one time step.
 *
 * Calculates the advection of the height field and the momentum update for the velocity field, storing the resulting derivatives in the provided arrays.
 *
 * @param t Current simulation time.
 * @param f Array of state variables (velocity components and height).
 * @param df Array to store computed time derivatives of the state variables.
 */
trace
void advance (double t, scalar * f, scalar * df)
{
  vector u = {f[0], f[1]}, du = {df[0], df[1]};
  scalar h = f[2], dh = df[2];

  advection_centered (h, u, dh);
  momentum (u, h, du);
}

void update (double t, scalar * f)
{
}

event defaults (i = 0)
{
  foreach()
    h[] = 1.;
}

event init (i = 0)
{
}

void run (void)
{
  init_grid (N);

  timer start = timer_start();
  iter = 0, t = 0;
  while (events (true)) {
    double dt = dtnext (timestep ());
#if 1
    advection_centered (h, u, hn);
    foreach()
      h[] += hn[]*dt;
    momentum (u, h, un);
    foreach_face()
      u.x[] += un.x[]*dt;
#else /* unstable! */
    scalar f[3] = { u, v, h };
    scalar df[2][3] = {{ un,  vn,  hn },
		       { un1, vn1, hn1 }};
    runge_kutta (2, t, dt, 3, f, df, advance, update);
#endif
    iter = inext, t = tnext;
  }
  timer_print (start, iter, 0);

  free_grid ();
}

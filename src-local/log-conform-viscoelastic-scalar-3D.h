/** Title: log-conform-viscoelastic-3D.h
# Version: 1.0
# Main feature 1: A exists in across the domain and relaxes according to \lambda. The stress only acts according to G.
# Main feature 2: This is the 3D implementation of [log-conform-viscoelastic-scalar-2D.h](log-conform-viscoelastic-scalar-2D.h).

# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Updated: Oct 19, 2024

# change log: Oct 19, 2024 (v1.0)
- 3D implementation
- scalar implementation
- Not tested yet.

*/

/** The code is same as http://basilisk.fr/src/log-conform.h but 
- written with G-\lambda formulation. 
- It also fixes the bug where [\sigma_p] = 0 & [\sigma_s] = \gamma\kappa instead of [\sigma_s+\sigma_p] = \gamma\kappa.
*/ 

/**
 * # TODO: 
 * axi compatibility is not there. 
 * This is full 3D. To keep the code easy to read, we will not implement axi compatibility. 
*/

#if AXI
#error "axi compatibility is not there. To keep the code easy to read, we will not implement axi compatibility just yet."
#endif

#include "bcg.h"

/*
conformation tensor */
// diagonal elements
scalar A11[], A22[], A33[];
// off-diagonal elements
scalar A12[], A13[], A23[];
/*
stress tensor */
// diagonal elements
scalar T11[], T22[], T33[];
// off-diagonal elements
scalar T12[], T13[], T23[];

event defaults (i = 0) {
  if (is_constant (a.x))
    a = new face vector;

  /*
  initialize A and T
  */
  for (scalar s in {A11, A22, A33}) {
    foreach () {
      s[] = 1.;
    }
  }
  for (scalar s in {T11, T12, T13, T22, T23, T33, A12, A13, A23}) {
    foreach(){
      s[] = 0.;
    }
  }

  for (scalar s in {A11, A22, A33, T11, T22, T33, A12, A13, A23, T12, T13, T23}) {
    if (s.boundary[left] != periodic_bc) {
      s[left] = neumann(0);
	    s[right] = neumann(0);
    }
    if (s.boundary[top] != periodic_bc) {
      s[top] = neumann(0);
	    s[bottom] = neumann(0);
    }
#if dimension == 3
    if (s.boundary[front] != periodic_bc) {
      s[front] = neumann(0);
	    s[back] = neumann(0);
    }
#endif
  }
}

/**
## Useful functions in 2D

The first step is to implement a routine to calculate the eigenvalues
and eigenvectors of the conformation tensor $\mathbf{A}$.

These structs ressemble Basilisk vectors and tensors but are just
arrays not related to the grid. */

#if dimension == 2
typedef struct { double x, y;}   pseudo_v;
typedef struct { pseudo_v x, y;} pseudo_t;

static void diagonalization_2D (pseudo_v * Lambda, pseudo_t * R, pseudo_t * A)
{
  /**
  The eigenvalues are saved in vector $\Lambda$ computed from the
  trace and the determinant of the symmetric conformation tensor
  $\mathbf{A}$. */

  if (sq(A->x.y) < 1e-15) {
    R->x.x = R->y.y = 1.;
    R->y.x = R->x.y = 0.;
    Lambda->x = A->x.x; Lambda->y = A->y.y;
    return;
  }

  double T = A->x.x + A->y.y; // Trace of the tensor
  double D = A->x.x*A->y.y - sq(A->x.y); // Determinant

  /**
  The eigenvectors, $\mathbf{v}_i$ are saved by columns in tensor
  $\mathbf{R} = (\mathbf{v}_1|\mathbf{v}_2)$. */

  R->x.x = R->x.y = A->x.y;
  R->y.x = R->y.y = -A->x.x;
  double s = 1.;
  for (int i = 0; i < dimension; i++) {
    double * ev = (double *) Lambda;
    ev[i] = T/2 + s*sqrt(sq(T)/4. - D);
    s *= -1;
    double * Rx = (double *) &R->x;
    double * Ry = (double *) &R->y;
    Ry[i] += ev[i];
    double mod = sqrt(sq(Rx[i]) + sq(Ry[i]));
    Rx[i] /= mod;
    Ry[i] /= mod;
  }
}
#endif

/*
Now this is the 3D implementation.
*/
#if dimension == 3

#include "eigen_decomposition.h"

typedef struct { double x, y, z; }   pseudo_v3d;
typedef struct { pseudo_v3d x, y, z; } pseudo_t3d;

static void diagonalization_3D (pseudo_v3d * Lambda, pseudo_t3d * R, pseudo_t3d * A)
{
  // Check if the matrix is already diagonal
  if (sq(A->x.y) + sq(A->x.z) + sq(A->y.z) < 1e-15) {
    R->x.x = R->y.y = R->z.z = 1.;
    R->y.x = R->x.y = R->z.x = R->x.z = R->z.y = R->y.z = 0.;
    Lambda->x = A->x.x; Lambda->y = A->y.y; Lambda->z = A->z.z;
    return;
  }

  // Compute eigenvalues using the eigen_decomposition function
  double matrix[3][3] = {{A->x.x, A->x.y, A->x.z},
                         {A->x.y, A->y.y, A->y.z},
                         {A->x.z, A->y.z, A->z.z}};
  double eigenvectors[3][3];
  double eigenvalues[3];
  
  compute_eigensystem_symmetric_3x3(matrix, eigenvectors, eigenvalues);

  // Store eigenvalues and eigenvectors
  Lambda->x = eigenvalues[0];
  Lambda->y = eigenvalues[1];
  Lambda->z = eigenvalues[2];

  R->x.x = eigenvectors[0][0]; R->x.y = eigenvectors[0][1]; R->x.z = eigenvectors[0][2];
  R->y.x = eigenvectors[1][0]; R->y.y = eigenvectors[1][1]; R->y.z = eigenvectors[1][2];
  R->z.x = eigenvectors[2][0]; R->z.y = eigenvectors[2][1]; R->z.z = eigenvectors[2][2];
}
#endif

/**
The stress tensor depends on previous instants and has to be
integrated in time. In the log-conformation scheme the advection of
the stress tensor is circumvented, instead the conformation tensor,
$\mathbf{A}$ (or more precisely the related variable $\Psi$) is
advanced in time.

In what follows we will adopt a scheme similar to that of [Hao \& Pan
(2007)](#hao2007). We use a split scheme, solving successively

a) the upper convective term:
$$
\partial_t \Psi = 2 \mathbf{B} + (\Omega \cdot \Psi -\Psi \cdot \Omega)
$$
b) the advection term:
$$
\partial_t \Psi + \nabla \cdot (\Psi \mathbf{u}) = 0
$$
c) the model term (but set in terms of the conformation 
tensor $\mathbf{A}$). In an Oldroyd-B viscoelastic fluid, the model is
$$ 
\partial_t \mathbf{A} = -\frac{\mathbf{f}_r (\mathbf{A})}{\lambda}
$$

The implementation below assumes that the values of $\Psi$ and
$\conform_p$ are never needed simultaneously. This means that $\conform_p$ can
be used to store (temporarily) the values of $\Psi$ (i.e. $\Psi$ is
just an alias for $\conform_p$). */

#if dimension == 2
event tracer_advection(i++)
{
  scalar Psi11 = A11;
  scalar Psi12 = A12;
  scalar Psi22 = A22;

  /**
  ### Computation of $\Psi = \log \mathbf{A}$ and upper convective term */

  foreach() {
    /**
      We assume that the stress tensor $\mathbf{\tau}_p$ depends on the
      conformation tensor $\mathbf{A}$ as follows
      $$
      \mathbf{\tau}_p = G_p (\mathbf{A}) =
      G_p (\mathbf{A} - I)
      $$
    */

    pseudo_t A;

    A.x.x = A11[];
    A.x.y = A12[];
    A.y.y = A22[];

    /**
    The conformation tensor is diagonalized through the
    eigenvector tensor $\mathbf{R}$ and the eigenvalues diagonal
    tensor, $\Lambda$. */

    pseudo_v Lambda;
    pseudo_t R;
    diagonalization_2D (&Lambda, &R, &A);
    
    /**
    $\Psi = \log \mathbf{A}$ is easily obtained after diagonalization, 
    $\Psi = R \cdot \log(\Lambda) \cdot R^T$. */
    
    Psi12[] = R.x.x*R.y.x*log(Lambda.x) + R.y.y*R.x.y*log(Lambda.y);
    Psi11[] = sq(R.x.x)*log(Lambda.x) + sq(R.x.y)*log(Lambda.y);
    Psi22[] = sq(R.y.y)*log(Lambda.y) + sq(R.y.x)*log(Lambda.x);

    /**
    We now compute the upper convective term $2 \mathbf{B} +
    (\Omega \cdot \Psi -\Psi \cdot \Omega)$.

    The diagonalization will be applied to the velocity gradient
    $(\nabla u)^T$ to obtain the antisymmetric tensor $\Omega$ and
    the traceless, symmetric tensor, $\mathbf{B}$. If the conformation
    tensor is $\mathbf{I}$, $\Omega = 0$ and $\mathbf{B}= \mathbf{D}$.  */

    pseudo_t B;
    double OM = 0.;
    if (fabs(Lambda.x - Lambda.y) <= 1e-20) {
      B.x.y = (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(4.*Delta); 
      foreach_dimension() 
        B.x.x = (u.x[1,0] - u.x[-1,0])/(2.*Delta);
    } else {
      pseudo_t M;
      foreach_dimension() {
        M.x.x = (sq(R.x.x)*(u.x[1] - u.x[-1]) +
        sq(R.y.x)*(u.y[0,1] - u.y[0,-1]) +
        R.x.x*R.y.x*(u.x[0,1] - u.x[0,-1] + 
        u.y[1] - u.y[-1]))/(2.*Delta);
        M.x.y = (R.x.x*R.x.y*(u.x[1] - u.x[-1]) + 
        R.x.y*R.y.x*(u.y[1] - u.y[-1]) +
        R.x.x*R.y.y*(u.x[0,1] - u.x[0,-1]) +
        R.y.x*R.y.y*(u.y[0,1] - u.y[0,-1]))/(2.*Delta);
      }
      double omega = (Lambda.y*M.x.y + Lambda.x*M.y.x)/(Lambda.y - Lambda.x);
      OM = (R.x.x*R.y.y - R.x.y*R.y.x)*omega;

      B.x.y = M.x.x*R.x.x*R.y.x + M.y.y*R.y.y*R.x.y;
      foreach_dimension()
        B.x.x = M.x.x*sq(R.x.x)+M.y.y*sq(R.x.y);	
    }

    /**
    We now advance $\Psi$ in time, adding the upper convective
    contribution. */

    double s = -Psi12[];
    Psi12[] += dt * (2. * B.x.y + OM * (Psi22[] - Psi11[]));
    s *= -1;
    Psi11[] += dt * 2. * (B.x.x + s * OM);
    s *= -1;
    Psi22[] += dt * 2. * (B.y.y + s * OM);
  }

  /**
  ### Advection of $\Psi$
 
  We proceed with step (b), the advection of the log of the
  conformation tensor $\Psi$. */

  advection ({Psi11, Psi12, Psi22}, uf, dt);

  /**
  ### Convert back to \conform_p */

  foreach() {
    /**
    It is time to undo the log-conformation, again by
    diagonalization, to recover the conformation tensor $\mathbf{A}$
    and to perform step (c).*/

    pseudo_t A = {{Psi11[], Psi12[]}, {Psi12[], Psi22[]}}, R;
    pseudo_v Lambda;
    diagonalization_2D (&Lambda, &R, &A);
    Lambda.x = exp(Lambda.x), Lambda.y = exp(Lambda.y);
    
    A.x.y = R.x.x*R.y.x*Lambda.x + R.y.y*R.x.y*Lambda.y;
    foreach_dimension()
      A.x.x = sq(R.x.x)*Lambda.x + sq(R.x.y)*Lambda.y;

    /**
    We perform now step (c) by integrating 
    $\mathbf{A}_t = -\mathbf{f}_r (\mathbf{A})/\lambda$ to obtain
    $\mathbf{A}^{n+1}$. This step is analytic,
    $$
    \int_{t^n}^{t^{n+1}}\frac{d \mathbf{A}}{\mathbf{I}- \mathbf{A}} = 
    \frac{\Delta t}{\lambda}
    $$
    */

    double intFactor = lambda[] != 0. ? exp(-dt/lambda[]): 0.;
    
    A.x.y *= intFactor;
    foreach_dimension()
      A.x.x = (1. - intFactor) + A.x.x*intFactor;

    /**
      Then the Conformation tensor $\mathcal{A}_p^{n+1}$ is restored from
      $\mathbf{A}^{n+1}$.  */
    
    A12[] = A.x.y;
    T12[] = Gp[]*A.x.y;
    A11[] = A.x.x;
    T11[] = Gp[]*(A.x.x - 1.);
    A22[] = A.y.y;
    T22[] = Gp[]*(A.y.y - 1.);
  }
}
#endif

#if dimension == 3
event tracer_advection(i++)
{
  scalar Psi11 = A11, Psi12 = A12, Psi13 = A13,
         Psi22 = A22, Psi23 = A23, Psi33 = A33;

  foreach() {
    pseudo_t3d A, R;
    pseudo_v3d Lambda;

    A.x.x = A11[]; A.x.y = A12[]; A.x.z = A13[];
    A.y.x = A12[]; A.y.y = A22[]; A.y.z = A23[];
    A.z.x = A13[]; A.z.y = A23[]; A.z.z = A33[];

    diagonalization_3D (&Lambda, &R, &A);
    
    // Compute Psi = log(A)
    Psi11[] = R.x.x*R.x.x*log(Lambda.x) + R.x.y*R.x.y*log(Lambda.y) + R.x.z*R.x.z*log(Lambda.z);
    Psi22[] = R.y.x*R.y.x*log(Lambda.x) + R.y.y*R.y.y*log(Lambda.y) + R.y.z*R.y.z*log(Lambda.z);
    Psi33[] = R.z.x*R.z.x*log(Lambda.x) + R.z.y*R.z.y*log(Lambda.y) + R.z.z*R.z.z*log(Lambda.z);
    Psi12[] = R.x.x*R.y.x*log(Lambda.x) + R.x.y*R.y.y*log(Lambda.y) + R.x.z*R.y.z*log(Lambda.z);
    Psi13[] = R.x.x*R.z.x*log(Lambda.x) + R.x.y*R.z.y*log(Lambda.y) + R.x.z*R.z.z*log(Lambda.z);
    Psi23[] = R.y.x*R.z.x*log(Lambda.x) + R.y.y*R.z.y*log(Lambda.y) + R.y.z*R.z.z*log(Lambda.z);

    // Compute B and Omega tensors (3D version)
    pseudo_t3d B, M, Omega;
    double OM[3] = {0., 0., 0.};

    if (fabs(Lambda.x - Lambda.y) <= 1e-20 && fabs(Lambda.y - Lambda.z) <= 1e-20) {
      // If all eigenvalues are equal, simplify calculations
      B.x.y = (u.y[1,0,0] - u.y[-1,0,0] + u.x[0,1,0] - u.x[0,-1,0])/(4.*Delta);
      B.x.z = (u.z[1,0,0] - u.z[-1,0,0] + u.x[0,0,1] - u.x[0,0,-1])/(4.*Delta);
      B.y.z = (u.z[0,1,0] - u.z[0,-1,0] + u.y[0,0,1] - u.y[0,0,-1])/(4.*Delta);
      B.x.x = (u.x[1,0,0] - u.x[-1,0,0])/(2.*Delta);
      B.y.y = (u.y[0,1,0] - u.y[0,-1,0])/(2.*Delta);
      B.z.z = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
      
      // Set Omega to zero for equal eigenvalues
      foreach_dimension()
        foreach_dimension()
          Omega.x.y = 0.;
    } else {
      // General case
      foreach_dimension() {
        M.x.x = (sq(R.x.x)*(u.x[1,0,0] - u.x[-1,0,0]) +
                sq(R.y.x)*(u.y[1,0,0] - u.y[-1,0,0]) +
                sq(R.z.x)*(u.z[1,0,0] - u.z[-1,0,0]))/(2.*Delta);
        M.x.y = (R.x.x*R.x.y*(u.x[1,0,0] - u.x[-1,0,0]) +
                R.y.x*R.y.y*(u.y[1,0,0] - u.y[-1,0,0]) +
                R.z.x*R.z.y*(u.z[1,0,0] - u.z[-1,0,0]))/(2.*Delta);
        M.x.z = (R.x.x*R.x.z*(u.x[1,0,0] - u.x[-1,0,0]) +
                R.y.x*R.y.z*(u.y[1,0,0] - u.y[-1,0,0]) +
                R.z.x*R.z.z*(u.z[1,0,0] - u.z[-1,0,0]))/(2.*Delta);
      }

      // Compute full Omega tensor
      double omega_xy = (Lambda.y*M.x.y + Lambda.x*M.y.x)/(Lambda.y - Lambda.x);
      double omega_xz = (Lambda.z*M.x.z + Lambda.x*M.z.x)/(Lambda.z - Lambda.x);
      double omega_yz = (Lambda.z*M.y.z + Lambda.y*M.z.y)/(Lambda.z - Lambda.y);

      Omega.x.y =  omega_xy; Omega.x.z =  omega_xz;
      Omega.y.x = -omega_xy; Omega.y.z =  omega_yz;
      Omega.z.x = -omega_xz; Omega.z.y = -omega_yz;
      Omega.x.x = Omega.y.y = Omega.z.z = 0.;

      // Compute B
      foreach_dimension() {
        B.x.x = M.x.x;
        B.x.y = 0.5 * (M.x.y + M.y.x);
        B.x.z = 0.5 * (M.x.z + M.z.x);
      }
    }

    Update Psi with upper convective term
    double dt2 = 2. * dt;
    Psi11[] += dt2 * (B.x.x + Omega.x.y*Psi12[] - Omega.y.x*Psi12[] + Omega.x.z*Psi13[] - Omega.z.x*Psi13[]);
    Psi22[] += dt2 * (B.y.y + Omega.y.x*Psi12[] - Omega.x.y*Psi12[] + Omega.y.z*Psi23[] - Omega.z.y*Psi23[]);
    Psi33[] += dt2 * (B.z.z + Omega.z.x*Psi13[] - Omega.x.z*Psi13[] + Omega.z.y*Psi23[] - Omega.y.z*Psi23[]);
    Psi12[] += dt2 * (B.x.y + Omega.x.x*Psi12[] - Omega.x.y*Psi11[] + Omega.x.y*Psi22[] - Omega.y.y*Psi12[] + Omega.x.z*Psi23[] - Omega.z.y*Psi13[]);
    Psi13[] += dt2 * (B.x.z + Omega.x.x*Psi13[] - Omega.x.z*Psi11[] + Omega.x.y*Psi23[] - Omega.y.z*Psi12[] + Omega.x.z*Psi33[] - Omega.z.z*Psi13[]);
    Psi23[] += dt2 * (B.y.z + Omega.y.x*Psi13[] - Omega.x.z*Psi12[] + Omega.y.y*Psi23[] - Omega.y.z*Psi22[] + Omega.y.z*Psi33[] - Omega.z.z*Psi23[]);
  }

  // Advection of Psi
  // advection ({Psi11, Psi12, Psi13, Psi22, Psi23, Psi33}, uf, dt);
  advection ({A11, A12, A13, A22, A23, A33}, uf, dt);

  // Convert back to A and T
  // foreach() {
  //   pseudo_t3d A, R;
  //   pseudo_v3d Lambda;
    
  //   A.x.x = Psi11[]; A.x.y = Psi12[]; A.x.z = Psi13[];
  //   A.y.x = Psi12[]; A.y.y = Psi22[]; A.y.z = Psi23[];
  //   A.z.x = Psi13[]; A.z.y = Psi23[]; A.z.z = Psi33[];

  //   diagonalization_3D (&Lambda, &R, &A);
    
  //   Lambda.x = exp(Lambda.x); Lambda.y = exp(Lambda.y); Lambda.z = exp(Lambda.z);

  //   // Reconstruct A
  //   A11[] = R.x.x*R.x.x*Lambda.x + R.x.y*R.x.y*Lambda.y + R.x.z*R.x.z*Lambda.z;
  //   A22[] = R.y.x*R.y.x*Lambda.x + R.y.y*R.y.y*Lambda.y + R.y.z*R.y.z*Lambda.z;
  //   A33[] = R.z.x*R.z.x*Lambda.x + R.z.y*R.z.y*Lambda.y + R.z.z*R.z.z*Lambda.z;
  //   A12[] = R.x.x*R.y.x*Lambda.x + R.x.y*R.y.y*Lambda.y + R.x.z*R.y.z*Lambda.z;
  //   A13[] = R.x.x*R.z.x*Lambda.x + R.x.y*R.z.y*Lambda.y + R.x.z*R.z.z*Lambda.z;
  //   A23[] = R.y.x*R.z.x*Lambda.x + R.y.y*R.z.y*Lambda.y + R.y.z*R.z.z*Lambda.z;

  //   // Apply relaxation
  //   double intFactor = lambda[] != 0. ? exp(-dt/lambda[]): 0.;
    
  //   foreach_dimension()
  //     A11[] = 1. + (A11[] - 1.)*intFactor;
  //   A12[] *= intFactor;
  //   A13[] *= intFactor;
  //   A23[] *= intFactor;

  //   // Compute T
  //   T11[] = Gp[]*(A11[] - 1.);
  //   T22[] = Gp[]*(A22[] - 1.);
  //   T33[] = Gp[]*(A33[] - 1.);
  //   T12[] = Gp[]*A12[];
  //   T13[] = Gp[]*A13[];
  //   T23[] = Gp[]*A23[];
  // }
}
#endif

/**
### Divergence of the viscoelastic stress tensor

The viscoelastic stress tensor $\mathbf{\tau}_p$ is defined at cell centers
while the corresponding force (acceleration) will be defined at cell
faces. Two terms contribute to each component of the momentum
equation. For example the $x$-component in Cartesian coordinates has
the following terms: $\partial_x \mathbf{\tau}_{p_{xx}} + \partial_y
\mathbf{\tau}_{p_{xy}}$. The first term is easy to compute since it can be
calculated directly from center values of cells sharing the face. The
other one is harder. It will be computed from vertex values. The
vertex values are obtained by averaging centered values.  Note that as
a result of the vertex averaging cells `[]` and `[-1,0]` are not
involved in the computation of shear. */

event acceleration (i++)
{
  face vector av = a;

  // foreach_face(x){
  //   if (fm.x[] > 1e-20) {
      
  //     double shearX = (T12[0,1]*cm[0,1] + T12[-1,1]*cm[-1,1] - 
  //     T12[0,-1]*cm[0,-1] - T12[-1,-1]*cm[-1,-1])/4.;
      
  //     av.x[] += (shearX + cm[]*T11[] - cm[-1]*T11[-1])*
  //     alpha.x[]/(sq(fm.x[])*Delta);
    
  //   }
  // }

  // foreach_face(y){
  //   if (fm.y[] > 1e-20) {

  //     double shearY = (T12[1,0]*cm[1,0] + T12[1,-1]*cm[1,-1] - 
  //     T12[-1,0]*cm[-1,0] - T12[-1,-1]*cm[-1,-1])/4.;
      
  //     av.y[] += (shearY + cm[]*T22[] - cm[0,-1]*T22[0,-1])*
  //     alpha.y[]/(sq(fm.y[])*Delta);

  //   }
  // }

}


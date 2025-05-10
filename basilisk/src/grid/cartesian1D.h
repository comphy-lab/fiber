typedef double real;

#define GRIDNAME "Cartesian 1D"
#define dimension 1
#define GHOSTS 1

#define _I     (point.i - 1)
#define _DELTA (1./point.n)

typedef struct {
  Grid g;
  char * d;
  int n;
} Cartesian;

struct _Point {
  int i, j, level, n;
};
static Point last_point;

#define cartesian ((Cartesian *)grid)

@define data(k,l,m) ((double *)&cartesian->d[(point.i + k)*datasize])
@define allocated(...) true

macro POINT_VARIABLES (Point point = point) { VARIABLES(); }

macro2 foreach (char flags = 0, Reduce reductions = None)
{
  OMP_PARALLEL (reductions) {
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point;
    point.n = cartesian->n;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 1; _k <= point.n; _k++) {
	point.i = _k;
	{...}
      }
  }
}
  
macro2 foreach_face_generic (char flags = 0, Reduce reductions = None,
				const char * order = "xyz")
{
  OMP_PARALLEL (reductions) {
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point;
    point.n = cartesian->n;
    int _k;
    OMP(omp for schedule(static))
      for (_k = 1; _k <= point.n + 1; _k++) {
	point.i = _k;
	{...}
      }
  }
}

/**
 * @brief Marks the current face as an x-direction face in 1D grid iterations.
 *
 * Intended for use within face iteration macros to indicate operations on x-direction faces.
 */
macro1 is_face_x() {{ int ig = -1; NOT_UNUSED(ig); {...} }}

// ghost cell coordinates for each direction
static int _ig[] = {1,-1};

/**
 * @brief Applies normal (non-periodic) box boundary conditions to ghost cells at a specified boundary.
 *
 * For each scalar in the provided list, sets the value in the ghost cell adjacent to the specified boundary using the scalar's boundary condition function. Only applies to scalars with non-periodic boundary conditions.
 *
 * @param b Pointer to the boundary structure specifying which boundary to apply.
 * @param list List of scalars to which the boundary condition should be applied.
 * @param l Grid refinement level (unused in this function).
 */

/**
 * @brief Applies box boundary conditions to scalars at a specified boundary.
 *
 * Separates scalars into centered and face-centered groups, applies their respective boundary conditions to ghost cells, and delegates normal face-centered scalars to the normal boundary handler. Ignores scalars with periodic boundary conditions.
 *
 * @param b Pointer to the boundary structure specifying which boundary to apply.
 * @param list List of scalars to which the boundary condition should be applied.
 * @param l Grid refinement level (unused in this function).
 */

/**
 * @brief Applies periodic boundary conditions in the x-direction for eligible scalars.
 *
 * For each scalar in the list with a periodic boundary condition at the right boundary, copies values from the opposite end of the grid to the ghost cells to enforce periodicity.
 *
 * @param b Pointer to the boundary structure (unused).
 * @param list List of scalars to which the periodic boundary should be applied.
 * @param l Grid refinement level (unused in this function).
 */

/**
 * @brief Frees all memory associated with the current grid and its boundaries.
 *
 * After calling this function, the global grid pointer is set to NULL.
 */

/**
 * @brief Sets all values of the specified scalars to a given value across the entire grid, including ghost cells.
 *
 * @param alist List of scalars whose values will be set.
 * @param val Value to assign to each scalar element.
 */

/**
 * @brief Initializes a 1D Cartesian grid of the specified size.
 *
 * Allocates memory for the grid and its data, initializes all scalar values to zero, and sets up box and periodic boundary conditions. If a grid of the same size already exists, no action is taken.
 *
 * @param n Number of grid points (excluding ghost cells).
 */

/**
 * @brief Reallocates the scalar data array to accommodate additional scalar fields.
 *
 * Adjusts the memory layout to insert new scalar data and updates the global data size accordingly.
 *
 * @param size Size in bytes to add for each new scalar.
 */

/**
 * @brief Locates the grid point corresponding to a given physical x-coordinate.
 *
 * Computes the grid index and level for the provided x-coordinate. The y and z coordinates are ignored in 1D.
 *
 * @param xp Physical x-coordinate.
 * @param yp Physical y-coordinate (ignored).
 * @param zp Physical z-coordinate (ignored).
 * @return Point structure with computed grid index and level.
 */

static void box_boundary_level_normal (const Boundary * b, scalar * list, int l)
{
  if (!list)
    return;
  
  int d = ((BoxBoundary *)b)->d;

  Point point;
  point.n = cartesian->n;
  int ig = _ig[d];
  assert (d <= left);
  point.i = d == right ? point.n + GHOSTS : GHOSTS;
  Point neighbor = {point.i + ig};
  for (scalar s in list) {
    scalar b = s.v.x;
    val(s,ig) = b.boundary[d] (point, neighbor, s, NULL);
  }
}

static double periodic_bc (Point point, Point neighbor, scalar s, bool * data);

static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  int d = ((BoxBoundary *)b)->d;
  scalar * centered = NULL, * normal = NULL;

  int component = d/2;
  for (scalar s in list)
    if (!is_constant(s) && s.boundary[d] != periodic_bc) {
      if (s.face) {
	if ((&s.d.x)[component]) {
	  scalar b = s.v.x;
	  if (b.boundary[d])
	    normal = list_add (normal, s);
	}
      }	
      else if (s.boundary[d])
	centered = list_add (centered, s);
    }

  if (centered) {
    Point point;
    point.n = cartesian->n;
    int ig = _ig[d];
    point.i = d == right ? point.n + GHOSTS - 1 : GHOSTS;
    Point neighbor = {point.i + ig};
    for (scalar s in centered)
      val(s,ig) = s.boundary[d] (point, neighbor, s, NULL);
    free (centered);
  }
    
  box_boundary_level_normal (b, normal, l);
  free (normal);
}

// periodic boundaries

static void periodic_boundary_level_x (const Boundary * b, scalar * list, int l)
{
  scalar * list1 = NULL;
  for (scalar s in list)
    if (!is_constant(s) && s.boundary[right] == periodic_bc)
      list1 = list_add (list1, s);
  if (!list1)
    return;

  Point point = *((Point *)grid);
  point.i = 0, point.n = N;
  for (int i = 0; i < GHOSTS; i++)
    for (scalar s in list1)
      s[i] = s[i + point.n];
  for (int i = point.n + GHOSTS; i < point.n + 2*GHOSTS; i++)
    for (scalar s in list1)
      s[i] = s[i - point.n];

  free (list1);
}

void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  free (cartesian->d);
  free (cartesian);
  grid = NULL;
}

@if TRASH
@ undef trash
@ define trash(list) reset(list, undefined)
@endif

void reset (void * alist, double val)
{
  scalar * list = (scalar *) alist;
  char * data = cartesian->d;
  for (int i = 0; i < cartesian->n + 2; i++, data += datasize) {
    double * v = (double *) data;
    for (scalar s in list)
      if (!is_constant(s))
	v[s.i] = val;
  }
}

void init_grid (int n)
{
  if (cartesian && n == cartesian->n)
    return;
  free_grid();
  Cartesian * p = qmalloc (1, Cartesian);
  size_t len = (n + 2)*datasize;
  p->n = N = n;
  p->d = qmalloc (len, char);
  /* trash the data just to make sure it's either explicitly
     initialised or never touched */
  double * v = (double *) p->d;
  for (int i = 0; i < len/sizeof(double); i++)
    v[i] = undefined;
  grid = (Grid *) p;
  reset (all, 0.);
  // box boundaries
  for (int d = 0; d < 2; d++) {
    BoxBoundary * box = qcalloc (1, BoxBoundary);
    box->d = d;
    Boundary * b = (Boundary *) box;
    b->level   = box_boundary_level;
    add_boundary (b);
  }
  // periodic boundaries
  Boundary * b = qcalloc (1, Boundary);
  b->level = periodic_boundary_level_x;
  add_boundary (b);
  // mesh size
  grid->n = grid->tn = n;
}

/**
 * @brief Expands the data array to accommodate additional scalar fields.
 *
 * Increases the storage size for each grid point by the specified number of bytes, shifting existing data as needed to make room for new scalar variables.
 *
 * @param size Number of bytes to add per grid point for new scalar data.
 */
void realloc_scalar (int size)
{
  Cartesian * p = cartesian;
  size_t len = (p->n + 2);
  qrealloc (p->d, len*(datasize + size), char);
  char * data = p->d + (len - 1)*datasize;
  for (int i = p->n + 1; i > 0; i--, data -= datasize)
    memmove (data + i*size, data, datasize);
  datasize += size;
}

Point locate (double xp = 0, double yp = 0, double zp = 0)
{
  Point point;
  point.n = cartesian->n;
  double a = (xp - X0)/L0*point.n;
  point.i = a + 1;
  point.level = (a > -0.5 && a < point.n + 0.5) ? 0 : - 1;
  return point;
}

#include "variables.h"
#include "cartesian-common.h"

macro2 foreach_vertex (char flags = 0, Reduce reductions = None)
{
  foreach_face_generic (flags, reductions) {
    int ig = -1; NOT_UNUSED(ig);
    {...}
  }
}

/**
 * @brief Initializes method pointers for the 1D Cartesian grid.
 *
 * Sets up the function table for grid operations specific to the 1D Cartesian implementation.
 */
void cartesian1D_methods()
{
  cartesian_methods();
}

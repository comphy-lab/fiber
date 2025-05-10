#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#if NVTX
 @include <nvtx3/nvToolsExt.h>
#endif

@define unmap(x,y)
@define trash(x)  // data trashing is disabled by default. Turn it on with
                  // -DTRASH=1

/**
 * @brief Placeholder macro for interpreter-overloaded constructs.
 *
 * This macro is intended to be replaced or expanded by the interpreter and has no effect in standard C compilation.
 */
auto macro2 BEGIN_FOREACH() {{...}}

#if _OPENMP
@ include <omp.h>
@ define OMP(x) Pragma(#x)

/**
 * @brief Temporarily disables OpenMP parallelization by redefining OMP macros as no-ops.
 *
 * This macro block undefines and redefines the OMP macro to disable OpenMP pragmas within its scope, then restores the original OMP macro definition after the block.
 */
macro OMP_SERIAL()
{
  @undef OMP
  @define OMP(x)
  ; // necessary so that the preproc above is included
  {...}
  @undef OMP
  @define OMP(x) Pragma(#x)
  ; // necessary so that the preproc above is included
}
   
#elif _MPI

@ /**
 * @brief Macro for OpenMP pragma insertion.
 *
 * Expands to an OpenMP pragma if OpenMP is enabled; otherwise, expands to nothing.
 */
define OMP(x)
macro OMP_SERIAL() {{...}}

@ include <mpi.h>
static int mpi_rank, mpi_npe;
@ define tid() mpi_rank
@ define pid() mpi_rank
@ define npe() mpi_npe

#else // not MPI, not OpenMP

@ /**
 * @brief Macro for OpenMP pragma insertion.
 *
 * Expands to the specified OpenMP pragma if OpenMP is enabled; otherwise, expands to nothing.
 *
 * @param x The OpenMP pragma directive (without the `#pragma` keyword).
 */
define OMP(x)
macro OMP_SERIAL() {{...}}

#endif // _MPI

#if _CADNA
# include <cadna.h>
#endif // CADNA

#if __cplusplus
# define delete delete_qcc
# define right right_qcc
# define left left_qcc
# define norm norm_qcc
# define new new_qcc
#endif // _cplusplus

@define _NVARMAX 65536
@define is_constant(v) ((v).i >= _NVARMAX)
@define constant(v) (is_constant(v) ? _constant[(v).i - _NVARMAX] : nodata)

@define systderr  stderr
@define systdout  stdout
#if _MPI
FILE * qstderr (void);
FILE * qstdout (void);
FILE * ferr = NULL, * fout = NULL;
@ /**
 * @brief Exits the program if running with multiple MPI processes.
 *
 * Prints an error message and terminates execution if the number of processing elements is greater than one, indicating the function is not yet compatible with MPI.
 */
def not_mpi_compatible()
do {
  if (npe() > 1) {
    fprintf (ferr, "%s() is not compatible with MPI (yet)\n", __func__);
    exit (1);
  }
} while(0)
@
@ define system(command) (pid() == 0 ? system(command) : 0)
#else
@ define qstderr() stderr
@ define qstdout() stdout
@ define ferr      stderr
@ define fout      stdout
@ define not_mpi_compatible()
#endif

/**
 * @brief Handles failed assertions by printing an error message and aborting execution.
 *
 * Prints the file name, line number, and failed condition to standard error, then terminates the program.
 *
 * @param file Source file where the assertion failed.
 * @param line Line number of the failed assertion.
 * @param cond String representation of the failed condition.
 */
static inline void qassert (const char * file, int line, const char * cond) {
  fprintf (stderr, "%s:%d: Assertion `%s' failed.\n", file, line, cond);
  abort();
}
#undef assert
#ifndef LINENO
# define LINENO __LINE__
#endif
#define assert(a) if (!(a)) qassert (__FILE__, LINENO, #a)

// Memory tracing
@define sysmalloc malloc
@define syscalloc calloc
@define sysrealloc realloc
@define sysfree free
@define systrdup strdup

#if MTRACE
# include "mtrace.h"
#else // !MTRACE
@ define pmalloc(s,func,file,line)    malloc(s)
@ define pcalloc(n,s,func,file,line)  calloc(n,s)
@ define prealloc(p,s,func,file,line) realloc(p,s)
@ define pfree(p,func,file,line)      free(p)
@ define pstrdup(s,func,file,line)    strdup(s)
#endif // !MTRACE

#define qrealloc(p, size, type) p = (type *) realloc (p, (size)*sizeof(type))
#define qmalloc(size, type) ((type *) malloc ((size)*sizeof(type)))
#define qcalloc(size, type) ((type *) calloc (size, sizeof(type)))

#include "array.h"

// Function tracing

#if TRACE == 1 // with Extrae library
#include <extrae_user_events.h>

typedef struct {
  Array index, stack;
  extrae_type_t type;
} Trace;

Trace trace_func     = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000010,
};

Trace trace_mpi_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000011,
};

/**
 * @brief Finds the index of a function name in an array or appends it if not found.
 *
 * Searches for the given function name in the provided array of strings. If found, returns its 1-based index. If not found, appends the function name to the array and returns the new length of the array.
 *
 * @param a Pointer to an Array structure containing function names as strings.
 * @param func The function name to search for or append.
 * @return int The 1-based index of the function name if found, or the new length of the array after appending.
 */
static int lookup_func (Array * a, const char * func)
{
  for (int i = 0; i < a->len/sizeof(char *); i++) {
    char * s = ((char **)a->p)[i];
    if (!strcmp (func, s))
      return i + 1;
  }
  char * s = strdup (func);
  array_append (a, &s, sizeof(char *));
  return a->len;
}

/**
 * @brief Records the entry of a function call in the tracing stack.
 *
 * Looks up the function identifier, emits a tracing event, and pushes the function onto the trace stack.
 *
 * @param t Pointer to the Trace structure managing the current tracing context.
 * @param func Name of the function being entered.
 */
static void trace_push (Trace * t, const char * func)
{
  int value = lookup_func (&t->index, func);
  Extrae_eventandcounters (t->type, value);
  array_append (&t->stack, &value, sizeof(int));
}

/**
 * @brief Pops the most recent function from the trace stack and updates the tracing event.
 *
 * Removes the top entry from the trace stack and signals the tracing system to update the current event context.
 *
 * @param t Pointer to the Trace structure managing the stack and event type.
 * @param func Name of the function being exited (unused in this implementation).
 */
static void trace_pop (Trace * t, const char * func)
{
  assert (t->stack.len > 0);
  t->stack.len -= sizeof(int);
  int value = t->stack.len > 0 ?
    ((int *)t->stack.p)[t->stack.len/sizeof(int) - 1] : 0;
  Extrae_eventandcounters (t->type, value);
}

/**
 * @brief Defines a new Extrae event type for function tracing.
 *
 * Registers a set of function names and their corresponding values as a new event type in the Extrae tracing system, using the provided description. The first event is labeled "OTHER" with value 0, followed by each function in the trace index.
 *
 * @param t Pointer to the Trace structure containing the function index.
 * @param description Description for the event type.
 */
static void trace_define (Trace * t, char * description)
{
  if (t->index.len > 0) {
    extrae_value_t values[t->index.len/sizeof(char *) + 1];
    char * names[t->index.len/sizeof(char *) + 1],
      ** func = (char **) t->index.p;
    names[0] = "OTHER";
    values[0] = 0;
    unsigned len = 1;
    for (int i = 0; i < t->index.len/sizeof(char *); i++, func++) {
      names[len] = *func;
      values[len++] = i + 1;
    }
    Extrae_define_event_type (&t->type, description, &len, values, names);
  }
}

/**
 * @brief Frees all memory associated with a Trace structure.
 *
 * Releases memory for function name strings, the function index array, and the call stack within the given Trace object.
 *
 * @param t Pointer to the Trace structure to be deallocated.
 */
static void trace_free (Trace * t)
{
  char ** func = (char **) t->index.p;
  for (int i = 0; i < t->index.len/sizeof(char *); i++, func++)
    free (*func);
  free (t->index.p);
  free (t->stack.p);
}

/**
 * @brief Finalizes and frees resources used for function tracing.
 *
 * This function defines and then releases tracing resources for both general and MPI-related Basilisk functions.
 */
static void trace_off()
{
  trace_define (&trace_func, "Basilisk functions");
  trace_define (&trace_mpi_func, "Basilisk functions (MPI-related)");
  trace_free (&trace_func);
  trace_free (&trace_mpi_func);
}
#if 0
#define TRACE_TYPE(func) (strncmp (func, "mpi_", 4) ?		\
			  &trace_func : &trace_func)
#else
#define TRACE_TYPE(func) &trace_func
#endif
@  define tracing(func, file, line)     trace_push (TRACE_TYPE(func), func)
@  define end_tracing(func, file, line) trace_pop (TRACE_TYPE(func), func)

#elif TRACE // built-in function tracing

typedef struct {
  char * func, * file;
  int line, calls;
  double total, self;
#if _MPI
  double min, max;
#endif // _MPI
} TraceIndex;
				      
struct {
  Array stack, index;
  double t0;
} Trace = {
  {NULL, 0, 0}, {NULL, 0, 0},
  -1
};

/**
 * @brief Records or updates profiling statistics for a function call in the tracing system.
 *
 * If the function, file, and line combination is new, adds a new entry to the trace index; otherwise, updates the existing entry with additional call count and timing information.
 *
 * @param func Name of the function being traced.
 * @param file Source file where the function is located.
 * @param line Line number in the source file.
 * @param total Total time spent in the function (including callees).
 * @param self Time spent exclusively in the function (excluding callees).
 */
static void trace_add (const char * func, const char * file, int line,
		       double total, double self)
{
  TraceIndex * t = (TraceIndex *) Trace.index.p;
  int i, len = Trace.index.len/sizeof(TraceIndex);
  for (i = 0; i < len; i++, t++)
    if (t->line == line && !strcmp (func, t->func) && !strcmp (file, t->file))
      break;
  if (i == len) {
    TraceIndex t = {strdup(func), strdup(file), line, 1, total, self};
    array_append (&Trace.index, &t, sizeof(TraceIndex));
  }
  else
    t->calls++, t->total += total, t->self += self;
}

/**
 * @brief Records the entry of a function for built-in tracing and profiling.
 *
 * Captures the current wall-clock time and pushes it onto the trace stack to mark the start of a function call. Also initiates an NVTX profiling range if enabled.
 *
 * @param func Name of the function being entered.
 * @param file Source file name where the function is called.
 * @param line Line number in the source file.
 */
static void tracing (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  if (Trace.t0 < 0)
    Trace.t0 = tv.tv_sec + tv.tv_usec/1e6;
  double t[2] = {(tv.tv_sec - Trace.t0) + tv.tv_usec/1e6, 0.};
  array_append (&Trace.stack, t, 2*sizeof(double));
#if 0
  fprintf (stderr, "trace %s:%s:%d t: %g sum: %g\n",
	   func, file, line, t[0], t[1]);
#endif
#if NVTX
  nvtxRangePush (func);
#endif
}

/**
 * @brief Marks the end of a traced function call and updates profiling statistics.
 *
 * Records the elapsed time since the corresponding tracing entry, updates the trace stack, and aggregates timing data for profiling. Also handles NVTX range pop if enabled.
 *
 * @param func Name of the function ending tracing.
 * @param file Source file name where tracing ends.
 * @param line Line number in the source file where tracing ends.
 */
static void end_tracing (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  double te = (tv.tv_sec - Trace.t0) + tv.tv_usec/1e6;
  double * t = (double *) Trace.stack.p;
  assert (Trace.stack.len >= 2*sizeof(double));
  t += Trace.stack.len/sizeof(double) - 2;
  Trace.stack.len -= 2*sizeof(double);
  double dt = te - t[0];
#if 0
  fprintf (stderr, "end trace %s:%s:%d ts: %g te: %g dt: %g sum: %g\n",
	   func, file, line, t[0], te, dt, t[1]);
#endif
  trace_add (func, file, line, dt, dt - t[1]);
  if (Trace.stack.len >= 2*sizeof(double)) {
    t -= 2;
    t[1] += dt;
  }
#if NVTX
  nvtxRangePop();
#endif
}

/**
 * @brief Comparison function for sorting TraceIndex structures by self time.
 *
 * Returns a value indicating the ordering of two TraceIndex pointers based on their self time.
 *
 * @param p1 Pointer to the first TraceIndex.
 * @param p2 Pointer to the second TraceIndex.
 * @return int Negative if t1->self < t2->self, zero if equal, positive otherwise.
 */
static int compar_self (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  return t1->self < t2->self;
}

#if _MPI
/**
 * @brief Comparison function for sorting TraceIndex structures by line number and file name.
 *
 * Compares two TraceIndex structures first by their line number, and if equal, by their file name using strcmp.
 *
 * @param p1 Pointer to the first TraceIndex structure.
 * @param p2 Pointer to the second TraceIndex structure.
 * @return Negative value if p1 < p2, zero if equal, positive if p1 > p2.
 */
static int compar_func (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  if (t1->line != t2->line)
    return t1->line < t2->line;
  return strcmp (t1->file, t2->file);
}
#endif

/**
 * @brief Prints a summary of function call profiling data exceeding a threshold.
 *
 * Outputs a table of traced function calls, including call counts, total and self times, and percentage of total time, to the specified file stream. Only functions whose self time exceeds the given percentage threshold are included. In MPI environments, aggregates and displays min/max statistics across all ranks.
 *
 * @param fp Output file stream.
 * @param threshold Minimum percentage of total self time for a function to be reported.
 */
void trace_print (FILE * fp, double threshold)
{
  int i, len = Trace.index.len/sizeof(TraceIndex);
  double total = 0.;
  TraceIndex * t;
  Array * index = array_new();
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    array_append (index, t, sizeof(TraceIndex)), total += t->self;
#if _MPI
  qsort (index->p, len, sizeof(TraceIndex), compar_func);
  double tot[len], self[len], min[len], max[len];
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    tot[i] = t->total, self[i] = t->self;
  MPI_Reduce (self,  min, len, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce (self,  max, len, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? self : MPI_IN_PLACE,
	      self, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? tot : MPI_IN_PLACE,
	      tot, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  total = 0.;
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    t->total = tot[i]/npe(), t->self = self[i]/npe(),
      t->max = max[i], t->min = min[i], total += t->self;
#endif // _MPI
  qsort (index->p, len, sizeof(TraceIndex), compar_self);
  fprintf (fp, "   calls    total     self   %% total   function\n");
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    if (t->self*100./total > threshold) {
      fprintf (fp, "%8d   %6.2f   %6.2f     %4.1f%%",
	       t->calls, t->total, t->self, t->self*100./total);
#if _MPI
      fprintf (fp, " (%4.1f%% - %4.1f%%)", t->min*100./total, t->max*100./total);
#endif
      fprintf (fp, "   %s():%s:%d\n", t->func, t->file, t->line);
    }
  fflush (fp);
  array_free (index);
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    t->calls = t->total = t->self = 0.;
}

/**
 * @brief Finalizes tracing and frees all associated resources.
 *
 * Outputs the collected trace summary and releases memory used for trace indices and stack data.
 */
static void trace_off()
{
  trace_print (fout, 0.);

  int i, len = Trace.index.len/sizeof(TraceIndex);
  TraceIndex * t;
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    free (t->func), free (t->file);

  free (Trace.index.p);
  Trace.index.p = NULL;
  Trace.index.len = Trace.index.max = 0;
  
  free (Trace.stack.p);
  Trace.stack.p = NULL;
  Trace.stack.len = Trace.stack.max = 0;
}

#else // disable tracing
@  define tracing(...)
@  define end_tracing(...)
#endif

// OpenMP / MPI
  
#if _OPENMP

@define tid() omp_get_thread_num()
@define pid() 0
@define npe() omp_get_num_threads()
@define mpi_all_reduce(v,type,op)
@define mpi_all_reduce_array(v,type,op,elem)

#elif _MPI

static bool in_prof = false;
static double prof_start, _prof;
@def prof_start(name)
  assert (!in_prof); in_prof = true;
  prof_start = MPI_Wtime();
@
@def prof_stop()
  assert (in_prof); in_prof = false;
  _prof = MPI_Wtime();
  mpi_time += _prof - prof_start;
@

#if FAKE_MPI
@define mpi_all_reduce(v,type,op)
@define mpi_all_reduce_array(v,type,op,elem)
#else /**
 * @brief Performs an MPI all-reduce operation on the provided buffer.
 *
 * Calls MPI_Allreduce to combine values from all processes and distribute the result to all processes in the communicator.
 *
 * @param sendbuf Pointer to the input buffer.
 * @param recvbuf Pointer to the output buffer where the result is stored.
 * @param count Number of elements in the buffer.
 * @param datatype MPI data type of buffer elements.
 * @param op MPI operation to perform (e.g., MPI_SUM, MPI_MAX).
 * @param comm MPI communicator.
 * @return int Result of the MPI_Allreduce call (MPI_SUCCESS on success).
 */
trace
int mpi_all_reduce0 (void *sendbuf, void *recvbuf, int count,
		     MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  return MPI_Allreduce (sendbuf, recvbuf, count, datatype, op, comm);
}

@/**
 * @brief Performs an MPI all-reduce operation on a single variable with profiling.
 *
 * Collectively reduces the value of a single variable across all MPI processes using the specified operation and data type, storing the result back in the original variable. Profiling is performed for the duration of the operation.
 */
def mpi_all_reduce(v,type,op) {
  prof_start ("mpi_all_reduce");
  union { int a; float b; double c;} global;
  mpi_all_reduce0 (&(v), &global, 1, type, op, MPI_COMM_WORLD);
  memcpy (&(v), &global, sizeof (v));
  prof_stop();
}
@

/**
 * @brief Performs an MPI all-reduce operation on an array.
 *
 * Applies the specified MPI reduction operation to an array of elements across all processes in MPI_COMM_WORLD, storing the result back into the input array. Aborts if the data type is not supported.
 *
 * @param v Pointer to the array to be reduced. The result overwrites the input.
 * @param datatype MPI data type of the array elements (supports MPI_DOUBLE, MPI_INT, MPI_LONG, MPI_C_BOOL, MPI_UNSIGNED_CHAR).
 * @param op MPI reduction operation (e.g., MPI_SUM, MPI_MAX).
 * @param elem Number of elements in the array.
 */
trace
void mpi_all_reduce_array (void * v, MPI_Datatype datatype, MPI_Op op, int elem)
{
  size_t size;
  if (datatype == MPI_DOUBLE) size = sizeof (double);
  else if (datatype == MPI_INT) size = sizeof (int);
  else if (datatype == MPI_LONG) size = sizeof (long);
  else if (datatype == MPI_C_BOOL) size = sizeof (bool);
  else if (datatype == MPI_UNSIGNED_CHAR) size = sizeof (unsigned char);
  else {
    fprintf (stderr, "unknown reduction type\n");
    fflush (stderr);
    abort();
  }
  void * global = malloc (elem*size), * tmp = malloc (elem*size);
  memcpy (tmp, v, elem*size);
  mpi_all_reduce0 (tmp, global, elem, datatype, op, MPI_COMM_WORLD);
  memcpy (v, global, elem*size);
  free (global), free (tmp);
}
#endif // !FAKE_MPI

@define QFILE /**
 * @brief Returns a FILE pointer for error output, redirected per MPI rank.
 *
 * On MPI rank 0, returns the system standard error stream. On other ranks, opens and returns a rank-specific log file for error output.
 *
 * @return FILE* Pointer to the appropriate error output stream for the current MPI rank.
 */
FILE // a dirty trick to avoid qcc 'static FILE *' rule

FILE * qstderr (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "log-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systderr;
  }
  return fp;
}

/**
 * @brief Returns a FILE pointer for standard output, redirecting to a rank-specific file in MPI environments.
 *
 * For MPI ranks greater than zero, output is written to a file named "out-<rank>". For rank zero or non-MPI environments, returns the standard output stream.
 *
 * @return FILE* Pointer to the appropriate output stream for the current process.
 */
FILE * qstdout (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "out-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systdout;
  }
  return fp;
}

/**
 * @brief Finalizes the MPI environment.
 *
 * Calls MPI_Finalize to cleanly shut down the MPI execution environment.
 */
static void finalize (void)
{
  MPI_Finalize();
}

/**
 * @brief Initializes the MPI environment and configures process-specific resources.
 *
 * Sets up MPI if not already initialized, assigns process rank and size, seeds the random number generator per rank, and redirects standard output and error streams for non-root ranks. Also configures memory tracing files per rank if enabled.
 */
void mpi_init()
{
  int initialized;
  MPI_Initialized (&initialized);
  if (!initialized) {
    MPI_Init (NULL, NULL);
    MPI_Comm_set_errhandler (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    atexit (finalize);
  }
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_npe);
  srand (mpi_rank + 1);
  if (ferr == NULL) {
    if (mpi_rank > 0) {
      ferr = fopen ("/dev/null", "w");
      fout = fopen ("/dev/null", "w");
    }
    else {
      ferr = systderr;
      fout = systdout;
    }
    char * etrace = getenv ("MALLOC_TRACE"), name[80];
    if (etrace && mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      setenv ("MALLOC_TRACE", name, 1);
    }
#if MTRACE == 1
    etrace = getenv ("MTRACE");
    if (!etrace)
      etrace = "mtrace";
    if (mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      pmtrace.fp = fopen (name, "w");
      pmtrace.fname = systrdup(name);
    }
    else {
      pmtrace.fp = fopen (etrace, "w");
      pmtrace.fname = systrdup(etrace);
    }
#endif
  }
}

#else // not MPI, not OpenMP

@define tid() 0
@define pid() 0
@define npe() 1
@define mpi_all_reduce(v,type,op)
@define mpi_all_reduce_array(v,type,op,elem)

#endif /**
 * @brief Placeholder macro for OpenMP parallel regions in non-parallel builds.
 *
 * Expands to an empty block when neither MPI nor OpenMP is enabled, allowing code to compile without parallelism.
 */

macro2 OMP_PARALLEL() {{...}}
@define OMP_PARALLEL(...) OMP(omp parallel S__VA_ARGS__)

@define NOT_UNUSED(x) (void)(x)

macro2 VARIABLES() { _CATCH; }
@define _index(a,m)    (a.i)
@define val(a,k,l,m)   data(k,l,m)[_index(a,m)]

double _val_higher_dimension = 0.;

/* undefined value */
/* Initialises unused memory with "signaling NaNs".  
 * This is probably not very portable, tested with
 * gcc (Debian 4.4.5-8) 4.4.5 on Linux 2.6.32-5-amd64.
 * This blog was useful:
 *   http://codingcastles.blogspot.co.nz/2008/12/nans-in-c.html 
 */
@if (_GNU_SOURCE || __APPLE__) && !_OPENMP && !_CADNA
double undefined;
@ if __APPLE__
@   include <stdint.h>
@   include "fp_osx.h"
@ endif
@if _GPU
@  /**
 * @brief Initializes the global undefined value as a signaling NaN and enables floating-point exceptions.
 *
 * Sets the global `undefined` variable to a signaling NaN bit pattern and enables floating-point exceptions for division by zero and invalid operations.
 */
define enable_fpe(flags)
@else
@  define enable_fpe(flags)  feenableexcept (flags)
@endif
@  define disable_fpe(flags) fedisableexcept (flags)
static void set_fpe (void) {
  int64_t lnan = 0x7ff0000000000001;
  assert (sizeof (int64_t) == sizeof (double));
  memcpy (&undefined, &lnan, sizeof (double));
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
@else // !((_GNU_SOURCE || __APPLE__) && !_OPENMP && !_CADNA && !_GPU)
@  /**
 * @brief No-op function for setting floating-point exception handling.
 *
 * This function is a placeholder and does not modify floating-point exception settings.
 */
define undefined ((double) DBL_MAX)
@  define enable_fpe(flags)
@  define disable_fpe(flags)
static void set_fpe (void) {}
@endif

// Pipes

static FILE ** qpopen_pipes = NULL;

/****
 * @brief Opens a pipe to a process, disabling pipes on non-root MPI ranks.
 *
 * On MPI rank 0, opens a pipe to the specified command using popen(). On other ranks, opens /dev/null instead. Tracks opened pipes for later cleanup.
 *
 * @param command Shell command to execute.
 * @param type Mode string as for popen(), e.g., "r" or "w".
 * @return FILE* Pointer to the opened pipe or file, or NULL on failure.
 */
FILE * qpopen (const char * command, const char * type)
{
  if (pid() > 0)
    return fopen ("/dev/null", type);
  FILE * fp = popen (command, type);
  if (fp) {
    FILE ** i = qpopen_pipes;
    int n = 0;
    while (i && *i) { n++; i++; }
    qrealloc (qpopen_pipes, n + 2, FILE *);
    qpopen_pipes[n] = fp;
    qpopen_pipes[n+1] = NULL;
  }
  return fp;
}

/**
 * @brief Closes a pipe opened by qpopen, tracking it for later cleanup.
 *
 * On non-root MPI ranks, closes the file as a regular file stream. On the root rank, marks the pipe as closed in the internal tracking list and closes it using pclose.
 *
 * @param fp Pointer to the FILE stream to close.
 * @return The result of fclose or pclose, depending on the process rank.
 */
int qpclose (FILE * fp)
{
  if (pid() > 0)
    return fclose (fp);
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i == fp)
      *i = (FILE *) 1;
    i++;
  }
  return pclose (fp);
}

/**
 * @brief Closes all pipes opened by qpopen and frees associated resources.
 *
 * Iterates through the list of open pipes created by qpopen, closes each one, and releases the memory used to track them.
 */
static void qpclose_all()
{
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i != (FILE *) 1)
      pclose (*i);
    i++;
  }
  free (qpopen_pipes);
  qpopen_pipes = NULL;
}

#define popen  qpopen
#define pclose qpclose

/**
 * @brief Opens a file with the process ID appended to its name.
 *
 * Constructs a filename by appending the current process ID to the given name, then opens the file with the specified mode.
 *
 * @param name Base filename to which the process ID will be appended.
 * @param mode File access mode (as in fopen).
 * @return FILE* Pointer to the opened file, or NULL if the file cannot be opened.
 */

FILE * lfopen (const char * name, const char * mode)
{
  char fname[80];
  sprintf (fname, "%s-%d", name, pid());
  return fopen (fname, mode);
}

#include "../ast/symbols.h"

enum typedef_kind_t {
  sym_SCALAR = sym_root + 1,
  sym_VECTOR,
  sym_TENSOR,
  sym_COORD,
  sym__COORD,
  sym_VEC4,
  sym_IVEC
};

@define attroffset(x) (offsetof(_Attributes,x))

/**
These are placeholders for internally-defined macros. */

typedef int Reduce;

macro2 foreach_face (char flags = 0, Reduce reductions = None,
			const char * order = "xyz")
{{...}}
/**
 * @brief Placeholder macro for Einstein summation notation.
 *
 * This macro is intended to be overloaded or implemented by an interpreter to support Einstein summation operations in mathematical expressions. By default, it has no effect.
 */
macro2 einstein_sum() {{...}}
/**
 * @brief Placeholder macro for diagonalization, intended for interpreter overloading.
 *
 * This macro is defined as a no-op in C and is expected to be overloaded by an interpreter
 * to provide diagonalization functionality when needed.
 */
macro2 diagonalize (int a) {{...}}

/**
Macros overloaded by the interpreter. */

@define dimensional(...)
#define show_dimension(...) show_dimension_internal (__VA_ARGS__ + 10293847566574839201.)
@define show_dimension_internal(...)
@define display_value(...)
@define interpreter_verbosity(...)

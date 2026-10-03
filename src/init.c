#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/RS.h>

/* Registered Fortran kernels (see src/panelio_smoothing.f90). Each entry
   declares the subroutine name, its entry point, the argument count and the
   per-argument SEXP types driving coercion in .Fortran (INTSXP = integer,
   REALSXP = double). */

extern void F77_SUB(smooth_kalman)(void *, void *, void *, void *,
                                   void *, void *, void *, void *);
extern void F77_SUB(smooth_kalman_nb)(void *, void *, void *, void *,
                                      void *, void *, void *, void *,
                                      void *, void *, void *, void *);
extern void F77_SUB(smooth_locreg)(void *, void *, void *, void *,
                                   void *, void *, void *, void *,
                                   void *, void *);
extern void F77_SUB(smooth_ewma)(void *, void *, void *, void *,
                                 void *, void *, void *, void *);
extern void F77_SUB(smooth_spline)(void *, void *, void *, void *,
                                   void *, void *, void *, void *);

/* demean2(mat, nobs, ntime, nfeat, maxiter, tol) */
extern void F77_SUB(demean2)(void *, void *, void *, void *,
                             void *, void *);
/* scm_effects(ctrl, tret, nobs, ncontrol, ntreated, nfeat, ntime, t0,
               weights, gap, att) */
extern void F77_SUB(scm_effects)(void *, void *, void *, void *,
                                 void *, void *, void *, void *,
                                 void *, void *, void *);
/* scm_simplex(ctrl, tret, nobs, ncontrol, ntreated, nfeat, ntime, t0, offset,
               maxiter, tol, weights, gap, att, wsum) */
extern void F77_SUB(scm_simplex)(void *, void *, void *, void *,
                                 void *, void *, void *, void *,
                                 void *, void *, void *, void *,
                                 void *, void *, void *);

/* smooth_kalman(x, nobs, nfeat, ntime, qratio, initvar, mode, y) */
static R_NativePrimitiveArgType kalman_types[] = {
    REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP
};

/* smooth_kalman_nb(x, nobs, nfeat, ntime, lib, phi, q, mode, n_em, dt,
                    use_dt, y) */
static R_NativePrimitiveArgType kalman_nb_types[] = {
    REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP,
    INTSXP, INTSXP, REALSXP, INTSXP, REALSXP
};

/* smooth_locreg(x, nobs, nfeat, ntime, kcode, degree, band, nbrk, brk, y) */
static R_NativePrimitiveArgType locreg_types[] = {
    REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP,
    INTSXP, INTSXP, REALSXP
};

/* smooth_ewma(x, nobs, nfeat, ntime, alpha, nbrk, brk, y) */
static R_NativePrimitiveArgType ewma_types[] = {
    REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP
};

/* smooth_spline(x, nobs, nfeat, ntime, lambda, nbrk, brk, y) */
static R_NativePrimitiveArgType spline_types[] = {
    REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP
};

/* demean2(mat, nobs, ntime, nfeat, maxiter, tol) */
static R_NativePrimitiveArgType demean2_types[] = {
    REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP
};

/* scm_effects(ctrl, tret, nobs, ncontrol, ntreated, nfeat, ntime, t0,
                weights, gap, att) */
static R_NativePrimitiveArgType scm_types[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,
    REALSXP, REALSXP, REALSXP
};

/* scm_simplex(ctrl, tret, nobs, ncontrol, ntreated, nfeat, ntime, t0, offset,
               maxiter, tol, weights, gap, att, wsum) */
static R_NativePrimitiveArgType scm_simplex_types[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP
};

static const R_FortranMethodDef FortranEntries[] = {
    {"smooth_kalman", (DL_FUNC) &F77_SUB(smooth_kalman), 8, kalman_types},
    {"smooth_kalman_nb", (DL_FUNC) &F77_SUB(smooth_kalman_nb), 12, kalman_nb_types},
    {"smooth_locreg", (DL_FUNC) &F77_SUB(smooth_locreg), 10, locreg_types},
    {"smooth_ewma",   (DL_FUNC) &F77_SUB(smooth_ewma),   8, ewma_types},
    {"smooth_spline", (DL_FUNC) &F77_SUB(smooth_spline), 8, spline_types},
    {"demean2",       (DL_FUNC) &F77_SUB(demean2),       6, demean2_types},
    {"scm_effects",   (DL_FUNC) &F77_SUB(scm_effects),  11, scm_types},
    {"scm_simplex",   (DL_FUNC) &F77_SUB(scm_simplex),  15, scm_simplex_types},
    {NULL, NULL, 0}
};

void R_init_panelio(DllInfo *dll) {
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

! =============================================================================
! src/panelio_causal.f90
! Causal-estimation kernels for the balanced abundance grid [observation, time,
! feature]. All per-feature iteration and the iterative two-way fixed-effect
! demeaning and the synthetic-control non-negative least squares live here in
! Fortran; the R layer dispatches one .Fortran call per operation, so there is
! no R-level loop.
!
!   demean2     iterative two-way (subject x time) within-demeaning in place
!   nnls_solve  Lawson-Hanson non-negative least squares (min ||A x - b||, x>=0)
!   scm_effects per-feature, per-treated-unit synthetic-control weights and
!               post-intervention gaps
! =============================================================================

subroutine cholsolve(a, n, b, x)
  ! Solve the symmetric positive-definite system A x = b by Cholesky.
  implicit none
  integer, intent(in)  :: n
  double precision, intent(in)  :: a(n, n), b(n)
  double precision, intent(out) :: x(n)
  double precision :: L(n, n), s
  integer :: i, j, k
  L = 0.0d0
  do j = 1, n
    s = a(j, j)
    do k = 1, j - 1
      s = s - L(j, k) ** 2
    end do
    L(j, j) = sqrt(max(s, 0.0d0))
    do i = j + 1, n
      s = a(i, j)
      do k = 1, j - 1
        s = s - L(i, k) * L(j, k)
      end do
      L(i, j) = s / L(j, j)
    end do
  end do
  x = 0.0d0
  do i = 1, n
    s = b(i)
    do k = 1, i - 1
      s = s - L(i, k) * x(k)
    end do
    x(i) = s / L(i, i)
  end do
  do i = n, 1, -1
    s = x(i)
    do k = i + 1, n
      s = s - L(k, i) * x(k)
    end do
    x(i) = s / L(i, i)
  end do
end subroutine cholsolve

subroutine nnls_solve(a, m, n, b, x)
  ! Lawson-Hanson NNLS: min ||A x - b||_2^2 subject to x >= 0.
  implicit none
  integer, intent(in)  :: m, n
  double precision, intent(in)  :: a(m, n), b(m)
  double precision, intent(out) :: x(n)
  double precision, parameter :: eps = 1.0d-12
  integer, parameter :: maxit = 4000
  logical :: passive(n)
  integer :: idx(n), np, jj
  double precision :: g(n), xp(n)
  double precision :: zeta, alpha, gmax
  integer :: i, j, k, it
  g = 0.0d0; xp = 0.0d0
  zeta = 0.0d0; alpha = 0.0d0; gmax = 0.0d0
  x = 0.0d0
  passive = .false.
  idx = 0; np = 0
  do it = 1, maxit
    ! dual gradient g = A^T (b - A x), computed only over the active set
    do j = 1, n
      if (passive(j)) cycle
      g(j) = 0.0d0
      do i = 1, m
        g(j) = g(j) + a(i, j) * b(i)
      end do
      do k = 1, n
        if (passive(k)) then
          do i = 1, m
            g(j) = g(j) - a(i, j) * a(i, k) * x(k)
          end do
        end if
      end do
    end do
    j = 0; gmax = 0.0d0
    do k = 1, n
      if (.not. passive(k)) then
        if (j == 0 .or. g(k) > gmax) then
          j = k; gmax = g(k)
        end if
      end if
    end do
    if (j == 0 .or. gmax <= eps) exit
    passive(j) = .true.
    np = 0
    do k = 1, n
      if (passive(k)) then
        np = np + 1
        idx(np) = k
      end if
    end do
    ! inner loop: solve on the passive set, move blocking variables out
    do
      call lsq_passive(xp)   ! xp(jj) = candidate for column idx(jj), jj = 1..np
      alpha = 1.0d30
      do jj = 1, np
        if (xp(jj) <= 0.0d0) then
          zeta = x(idx(jj)) / (x(idx(jj)) - xp(jj))
          if (zeta < alpha) alpha = zeta
        end if
      end do
      if (alpha > 1.0d29) then
        do jj = 1, np
          x(idx(jj)) = xp(jj)
        end do
        exit
      else
        do jj = 1, np
          x(idx(jj)) = x(idx(jj)) + alpha * (xp(jj) - x(idx(jj)))
        end do
        do jj = 1, np
          if (x(idx(jj)) <= eps) then
            passive(idx(jj)) = .false.
            x(idx(jj)) = 0.0d0
          end if
        end do
        np = 0
        do k = 1, n
          if (passive(k)) then
            np = np + 1
            idx(np) = k
          end if
        end do
      end if
    end do
  end do
contains
  subroutine lsq_passive(xp)
    ! Solve min ||A_P xp - b||^2 over the passive columns (given by the host
    ! idx(1..np)) by a numerically stable modified Gram-Schmidt QR, robust to
    ! rank deficiency. The returned xp is indexed by passive-slot position jj,
    ! i.e. xp(jj) corresponds to column idx(jj).
    double precision, intent(out) :: xp(n)
    integer :: nn, jj, i, k
    double precision :: ap(m, n), R(n, n), qb(n), s
    nn = np
    ap = 0.0d0; R = 0.0d0; qb = 0.0d0; s = 0.0d0
    if (nn == 0) then
      xp = 0.0d0
      return
    end if
    do jj = 1, nn
      do i = 1, m
        ap(i, jj) = a(i, idx(jj))
      end do
    end do
    R = 0.0d0
    do jj = 1, nn
      s = 0.0d0
      do i = 1, m
        s = s + ap(i, jj) * ap(i, jj)
      end do
      R(jj, jj) = sqrt(s)
      if (R(jj, jj) < 1.0d-11) then
        R(jj, jj) = 0.0d0      ! numerically zero column: leave it out
        qb(jj) = 0.0d0
        cycle
      end if
      do i = 1, m
        ap(i, jj) = ap(i, jj) / R(jj, jj)
      end do
      s = 0.0d0
      do i = 1, m
        s = s + ap(i, jj) * b(i)
      end do
      qb(jj) = s
      do k = jj + 1, nn
        s = 0.0d0
        do i = 1, m
          s = s + ap(i, jj) * ap(i, k)
        end do
        R(jj, k) = s
        do i = 1, m
          ap(i, k) = ap(i, k) - s * ap(i, jj)
        end do
      end do
    end do
    xp = 0.0d0
    do jj = nn, 1, -1
      if (R(jj, jj) == 0.0d0) cycle
      s = qb(jj)
      do k = jj + 1, nn
        s = s - R(jj, k) * xp(k)
      end do
      xp(jj) = s / R(jj, jj)
    end do
  end subroutine lsq_passive
end subroutine nnls_solve

subroutine demean2(mat, nobs, ntime, nfeat, maxiter, tol)
  ! Iterative two-way within-demeaning (subject and time fixed effects) of each
  ! feature grid in place. Each sweep subtracts time means then subject means
  ! until the largest removed mean is below tol or maxiter is reached.
  implicit none
  integer, intent(in)  :: nobs, ntime, nfeat, maxiter
  double precision, intent(inout) :: mat(nobs, ntime, nfeat)
  double precision, intent(in)  :: tol
  integer :: iter, i, t, f
  double precision :: m, mx
  do f = 1, nfeat
    do iter = 1, maxiter
      mx = 0.0d0
      do t = 1, ntime
        m = 0.0d0
        do i = 1, nobs
          m = m + mat(i, t, f)
        end do
        m = m / dble(nobs)
        if (abs(m) > mx) mx = abs(m)
        do i = 1, nobs
          mat(i, t, f) = mat(i, t, f) - m
        end do
      end do
      do i = 1, nobs
        m = 0.0d0
        do t = 1, ntime
          m = m + mat(i, t, f)
        end do
        m = m / dble(ntime)
        if (abs(m) > mx) mx = abs(m)
        do t = 1, ntime
          mat(i, t, f) = mat(i, t, f) - m
        end do
      end do
      if (mx < tol) exit
    end do
  end do
end subroutine demean2

subroutine scm_effects(ctrl, tret, nobs, ncontrol, ntreated, nfeat, ntime, &
                       t0, weights, gap, att)
  ! Synthetic control per feature and per treated unit. For each treated unit
  ! the non-negative weights over control units minimise the squared mismatch
  ! of the pre-intervention trajectory; the post-intervention gap is the
  ! treated trajectory minus the weighted control trajectory. att is the mean
  ! over treated units of the mean post-intervention gap per feature.
  implicit none
  integer, intent(in) :: nobs, ncontrol, ntreated, nfeat, ntime, t0
  double precision, intent(in)  :: ctrl(ncontrol, ntime, nfeat)
  double precision, intent(in)  :: tret(ntreated, ntime, nfeat)
  double precision, intent(out) :: weights(ntreated, ncontrol, nfeat)
  double precision, intent(out) :: gap(ntreated, ntime - t0 + 1, nfeat)
  double precision, intent(out) :: att(nfeat)
  integer :: f, i, c, tt, tpre, tpost
  ! Automatic arrays whose leading dimension is t0-1 (tpre) so that
  ! nnls_solve(a, tpre, ncontrol, ...) sees the correct memory layout.
  double precision :: a(t0 - 1, ncontrol), b(t0 - 1), x(ncontrol), s
  a = 0.0d0; b = 0.0d0; x = 0.0d0; s = 0.0d0
  tpre = t0 - 1
  tpost = ntime - t0 + 1
  att = 0.0d0
  do f = 1, nfeat
    do i = 1, ntreated
      do c = 1, ncontrol
        do tt = 1, tpre
          a(tt, c) = ctrl(c, tt, f)
        end do
      end do
      do tt = 1, tpre
        b(tt) = tret(i, tt, f)
      end do
      call nnls_solve(a, tpre, ncontrol, b, x)
      do c = 1, ncontrol
        weights(i, c, f) = x(c)
      end do
      do tt = 1, tpost
        s = tret(i, t0 + tt - 1, f)
        do c = 1, ncontrol
          s = s - ctrl(c, t0 + tt - 1, f) * x(c)
        end do
        gap(i, tt, f) = s
      end do
    end do
    s = 0.0d0
    do i = 1, ntreated
      do tt = 1, tpost
        s = s + gap(i, tt, f)
      end do
    end do
    att(f) = s / dble(ntreated * tpost)
  end do
end subroutine scm_effects


! ---------------------------------------------------------------------------
! Exactly solve min ||A w - b||^2 subject to w >= 0 and sum(w) = 1 by an
! active-set KKT method. At the optimum there is a multiplier lambda with
!     grad_i = A_i' (A w - b) = lambda      for active (w_i > 0)
!     grad_i                              <= lambda      for inactive (w_i = 0)
! The method starts with every donor active, solves the reduced KKT system
!   [ A_S' A_S   1 ] [w_S]   [ A_S' b ]
!   [   1'       0 ] [lambda] = [   1   ]
! for the active set S, drops a negative weight if one appears, and admits a
! violating inactive donor, until the KKT conditions hold to the tolerance.
! ---------------------------------------------------------------------------
subroutine solve_simplex(a, m, n, b, tol, w)
  implicit none
  integer, intent(in) :: m, n
  double precision, intent(in)  :: a(m, n), b(m), tol
  double precision, intent(out) :: w(n)
  integer :: nS, i, j, k, it, p
  double precision :: H(n, n), g(n)
  double precision :: lam, grad, mx, hm, a1, b1
  double precision :: v(n), u(n), onev(n)
  double precision, parameter :: tiny_d = 1.0d-14
  logical :: active(n)
  integer :: idx(n)
  ! H = A'A, g = A'b
  do i = 1, n
    g(i) = 0.0d0
    do k = 1, m
      g(i) = g(i) + a(k, i) * b(k)
    end do
    do j = 1, n
      H(i, j) = 0.0d0
      do k = 1, m
        H(i, j) = H(i, j) + a(k, i) * a(k, j)
      end do
    end do
  end do
  ! feasible build-up start: put all weight on the donor with the strongest
  ! positive alignment with the target (largest A'b), so the active set stays
  ! small and the reduced Hessian H_S remains non-singular even when the number
  ! of pre-intervention rows tpre is smaller than the donor-pool size ncontrol
  w = 0.0d0
  do i = 1, n
    onev(i) = 1.0d0
  end do
  ! start donor: strongest alignment with the target among non-degenerate
  ! columns (H(i,i) > tiny), ties broken toward the larger reduced diagonal so
  ! the single-donor active set is never singular. When the target pre-period is
  ! all-zero (a structural zero), g is identically zero for every donor; a naive
  ! max-g scan would then select the last column, which can be a time-constant
  ! donor with H(i,i) = 0 and collapse the weights to zero.
  mx = -1.0d30
  hm = -1.0d30
  p = 0
  do i = 1, n
    if (H(i, i) > tiny_d) then
      if (g(i) > mx .or. (abs(g(i) - mx) <= tiny_d .and. H(i, i) > hm)) then
        mx = g(i)
        hm = H(i, i)
        p = i
      end if
    end if
  end do
  if (p == 0) then
    ! every donor is time-constant on this feature: no matching signal exists,
    ! so fall back to uniform convex weights (sum(w) = 1) instead of all zeros.
    do i = 1, n
      w(i) = 1.0d0 / dble(n)
    end do
    return
  end if
  active = .false.
  active(p) = .true.
  nS = 1
  do it = 1, n * 4 + 2
    ! rebuild the active index list
    p = 0
    do i = 1, n
      if (active(i)) then
        p = p + 1
        idx(p) = i
      end if
    end do
    nS = p
    ! closed-form KKT: v = H_S^{-1} 1, u = H_S^{-1} g_S, lambda = (1'u - 1)/(1'v)
    call pdsolve(H, idx, nS, n, onev, v)
    call pdsolve(H, idx, nS, n, g, u)
    a1 = 0.0d0
    b1 = 0.0d0
    do i = 1, nS
      a1 = a1 + v(i)
      b1 = b1 + u(i)
    end do
    lam = (b1 - 1.0d0) / max(a1, tiny_d)
    do i = 1, n
      w(i) = 0.0d0
    end do
    do i = 1, nS
      w(idx(i)) = u(i) - lam * v(i)
    end do
    ! drop the most negative active weight if any
    mx = -tol
    p = 0
    do i = 1, nS
      if (w(idx(i)) < mx) then
        mx = w(idx(i))
        p = i
      end if
    end do
    if (p /= 0) then
      active(idx(p)) = .false.
      cycle
    end if
    ! admit the most violating inactive donor (grad_i > -lambda, above the
    ! common active gradient c = -lambda)
    mx = -1.0d30
    p = 0
    do i = 1, n
      if (.not. active(i)) then
        grad = 0.0d0
        do j = 1, n
          grad = grad + H(i, j) * w(j)
        end do
        grad = grad - g(i)
        if (grad > -lam + tol .and. grad > mx) then
          mx = grad
          p = i
        end if
      end if
    end do
    if (p == 0) exit   ! KKT satisfied
    active(p) = .true.
  end do
contains
  subroutine pdsolve(H, idx, ns, n, rhs, sol)
    ! Solve H_S x = rhs over the active columns H(idx,idx) by Cholesky, with
    ! zero/negative pivots skipped so that a numerically rank-deficient active
    ! set returns a finite least-norm-style solution rather than NaN.
    integer, intent(in) :: idx(n), ns, n
    double precision, intent(in) :: H(n, n), rhs(n)
    double precision, intent(out) :: sol(ns)
    double precision :: L(ns, ns), s, diag
    integer :: i, j, k
    L = 0.0d0
    do j = 1, ns
      s = H(idx(j), idx(j))
      do k = 1, j - 1
        s = s - L(j, k) * L(j, k)
      end do
      diag = sqrt(max(s, 0.0d0))
      if (diag < 1.0d-14) then
        L(j, j) = 0.0d0
      else
        L(j, j) = diag
      end if
      do i = j + 1, ns
        s = H(idx(i), idx(j))
        do k = 1, j - 1
          s = s - L(i, k) * L(j, k)
        end do
        if (L(j, j) > 0.0d0) then
          L(i, j) = s / L(j, j)
        else
          L(i, j) = 0.0d0
        end if
      end do
    end do
    sol = 0.0d0
    do i = 1, ns
      s = rhs(idx(i))
      do k = 1, i - 1
        s = s - L(i, k) * sol(k)
      end do
      if (L(i, i) > 0.0d0) sol(i) = s / L(i, i)
    end do
    do i = ns, 1, -1
      s = sol(i)
      do k = i + 1, ns
        s = s - L(k, i) * sol(k)
      end do
      if (L(i, i) > 0.0d0) sol(i) = s / L(i, i)
    end do
  end subroutine pdsolve
end subroutine solve_simplex

! ---------------------------------------------------------------------------
! Exactly-convex synthetic control (the paper's default): per feature and per
! treated unit, the weights over control units solve
!     min ||X_pre w + mu * 1 - Y_pre||^2  s.t.  w >= 0, sum(w) = 1
! by an exact active-set KKT solve over the simplex, with the scalar level offset mu
! (centred log-ratio translation) optionally absorbed by demeaning the pre
! equations. The counterfactual always lies in the convex hull of the donor
! trajectories; the post-intervention gap is Y_post - (X_post w + mu). wsum
! carries the exact weight sum (1) so feasibility is reportable.
! ---------------------------------------------------------------------------
subroutine scm_simplex(ctrl, tret, nobs, ncontrol, ntreated, nfeat, ntime, &
                       t0, offset, maxiter, tol, weights, gap, att, wsum)
  implicit none
  integer, intent(in) :: nobs, ncontrol, ntreated, nfeat, ntime, t0, offset
  integer, intent(in) :: maxiter
  double precision, intent(in)  :: tol
  double precision, intent(in)  :: ctrl(ncontrol, ntime, nfeat)
  double precision, intent(in)  :: tret(ntreated, ntime, nfeat)
  double precision, intent(out) :: weights(ntreated, ncontrol, nfeat)
  double precision, intent(out) :: gap(ntreated, ntime - t0 + 1, nfeat)
  double precision, intent(out) :: att(nfeat)
  double precision, intent(out) :: wsum(ntreated, nfeat)
  integer :: f, i, c, tt, tpre, tpost
  double precision :: a(t0 - 1, ncontrol), b(t0 - 1)
  double precision :: ac(t0 - 1, ncontrol), bc(t0 - 1)
  double precision :: am(ncontrol), w(ncontrol)
  double precision :: s, mu, bmean
  tpre = t0 - 1
  tpost = ntime - t0 + 1
  att = 0.0d0
  do f = 1, nfeat
    do i = 1, ntreated
      do c = 1, ncontrol
        do tt = 1, tpre
          a(tt, c) = ctrl(c, tt, f)
        end do
      end do
      do tt = 1, tpre
        b(tt) = tret(i, tt, f)
      end do
      ! Optional level offset: demean the rows to remove the free scalar mu.
      if (offset == 1) then
        bmean = 0.0d0
        do tt = 1, tpre
          bmean = bmean + b(tt)
        end do
        bmean = bmean / dble(tpre)
        do c = 1, ncontrol
          am(c) = 0.0d0
          do tt = 1, tpre
            am(c) = am(c) + a(tt, c)
          end do
          am(c) = am(c) / dble(tpre)
          do tt = 1, tpre
            ac(tt, c) = a(tt, c) - am(c)
          end do
        end do
        do tt = 1, tpre
          bc(tt) = b(tt) - bmean
        end do
      else
        bmean = 0.0d0
        am = 0.0d0
        do tt = 1, tpre
          bc(tt) = b(tt)
        end do
        do c = 1, ncontrol
          do tt = 1, tpre
            ac(tt, c) = a(tt, c)
          end do
        end do
      end if
      ! exactly-convex weights solved as a small QP over the simplex by an
      ! active-set KKT solver (w >= 0, sum(w) = 1), exact to the tolerance
      call solve_simplex(ac, tpre, ncontrol, bc, tol, w)
      wsum(i, f) = 0.0d0
      do c = 1, ncontrol
        wsum(i, f) = wsum(i, f) + w(c)
      end do
      do c = 1, ncontrol
        weights(i, c, f) = w(c)
      end do
      ! level offset for the counterfactual (0 unless the rows were demeaned)
      mu = bmean
      do c = 1, ncontrol
        mu = mu - am(c) * w(c)
      end do
      do tt = 1, tpost
        s = tret(i, t0 + tt - 1, f)
        do c = 1, ncontrol
          s = s - ctrl(c, t0 + tt - 1, f) * w(c)
        end do
        gap(i, tt, f) = s - mu
      end do
    end do
    s = 0.0d0
    do i = 1, ntreated
      do tt = 1, tpost
        s = s + gap(i, tt, f)
      end do
    end do
    att(f) = s / dble(ntreated * tpost)
  end do
end subroutine scm_simplex

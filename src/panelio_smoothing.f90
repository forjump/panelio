! =============================================================================
! src/panelio_smoothing.f90
! Time-axis smoothing kernels for the abundance array [observation, feature,
! time]. All iteration over time (and over the observation x feature grid) is
! done here in Fortran; the R layer dispatches exactly one .Fortran call per
! method, so there is no R-level loop. Every method except the Kalman family
! accepts breakpoints (1-based time indices that start a new segment): each
! segment is smoothed independently so an intervention effect at the breakpoint
! is never smoothed away. The Kalman family ignores breakpoints by design.
! =============================================================================

subroutine build_segments(ntime, nbrk, brk, nseg, lo, hi)
  ! Splits [1, ntime] at breakpoints brk(1) < brk(2) < ... < brk(nbrk), each of
  ! which is the 1-based index that STARTS a new (post-intervention) segment.
  ! Segments: [1, b1-1], [b1, b2-1], ..., [b_{nbrk-1}, b_nbrk-1], [b_nbrk, ntime].
  implicit none
  integer, intent(in)  :: ntime, nbrk
  integer, intent(in)  :: brk(*)
  integer, intent(out) :: nseg
  integer, intent(out) :: lo(*), hi(*)
  integer :: s
  nseg = nbrk + 1
  if (nbrk >= 1) then
    lo(1) = 1
    hi(1) = brk(1) - 1
    do s = 2, nbrk
      lo(s) = brk(s - 1)
      hi(s) = brk(s) - 1
    end do
    lo(nseg) = brk(nbrk)
    hi(nseg) = ntime
  else
    lo(1) = 1
    hi(1) = ntime
  end if
  ! enforce a consistent, non-overlapping ordering
  do s = 2, nseg
    if (lo(s) <= hi(s - 1)) lo(s) = hi(s - 1) + 1
  end do
  do s = 1, nseg
    if (hi(s) < lo(s)) hi(s) = lo(s) - 1   ! mark an empty segment
  end do
end subroutine build_segments

subroutine solve_poly(a, b, m, sol)
  ! Solve the small symmetric normal-equation system A sol = b (m <= 4) by
  ! Gaussian elimination with partial pivoting.
  implicit none
  integer, intent(in)  :: m
  double precision, intent(inout) :: a(m, m), b(m)
  double precision, intent(out)   :: sol(m)
  double precision :: am(m, m + 1), pivot, factor
  integer :: i, j, k, p
  do i = 1, m
    do j = 1, m
      am(i, j) = a(i, j)
    end do
    am(i, m + 1) = b(i)
  end do
  do k = 1, m
    p = k
    do i = k + 1, m
      if (abs(am(i, k)) > abs(am(p, k))) p = i
    end do
    if (p /= k) then
      do j = k, m + 1
        pivot = am(k, j); am(k, j) = am(p, j); am(p, j) = pivot
      end do
    end if
    do i = k + 1, m
      factor = am(i, k) / am(k, k)
      am(i, k) = 0.0d0
      do j = k + 1, m + 1
        am(i, j) = am(i, j) - factor * am(k, j)
      end do
    end do
  end do
  sol(m) = am(m, m + 1) / am(m, m)
  do i = m - 1, 1, -1
    sol(i) = am(i, m + 1)
    do j = i + 1, m
      sol(i) = sol(i) - am(i, j) * sol(j)
    end do
    sol(i) = sol(i) / am(i, i)
  end do
end subroutine solve_poly

! ---------------------------------------------------------------------------
! Kalman family: forward filter (EKF/KF) or RTS smoother. Breakpoints ignored.
! ---------------------------------------------------------------------------
subroutine smooth_kalman(x, nobs, nfeat, ntime, qratio, initvar, mode, y)
  implicit none
  integer, intent(in)  :: nobs, nfeat, ntime, mode
  double precision, intent(in)  :: x(nobs, nfeat, ntime)
  double precision, intent(in)  :: qratio, initvar
  double precision, intent(out) :: y(nobs, nfeat, ntime)
  double precision :: series(ntime)
  double precision :: xf(ntime), pf(ntime), ppred(ntime)
  double precision :: r, q, meanv, var, k, pp, xsm, psm, gain
  integer :: i, j, t
  do j = 1, nfeat
    do i = 1, nobs
      do t = 1, ntime
        series(t) = x(i, j, t)
      end do
      meanv = 0.0d0
      do t = 1, ntime
        meanv = meanv + series(t)
      end do
      meanv = meanv / dble(ntime)
      var = 0.0d0
      do t = 1, ntime
        var = var + (series(t) - meanv) ** 2
      end do
      if (ntime > 1) var = var / dble(ntime - 1)
      if (var <= 0.0d0) var = 1.0d0
      r = var
      q = qratio * var
      xf(1) = series(1)
      pf(1) = initvar * var
      y(i, j, 1) = xf(1)
      do t = 2, ntime
        ppred(t) = pf(t - 1) + q
        k = ppred(t) / (ppred(t) + r)
        xf(t) = xf(t - 1) + k * (series(t) - xf(t - 1))
        pf(t) = (1.0d0 - k) * ppred(t)
        y(i, j, t) = xf(t)
      end do
      if (mode == 1) then
        xsm = xf(ntime)
        psm = pf(ntime)
        y(i, j, ntime) = xsm
        do t = ntime - 1, 1, -1
          gain = pf(t) / ppred(t + 1)
          xsm = xf(t) + gain * (xsm - xf(t))
          psm = pf(t) + gain * (psm - ppred(t + 1)) * gain
          y(i, j, t) = xsm
        end do
      end if
    end do
  end do
end subroutine smooth_kalman

! ---------------------------------------------------------------------------
! Local polynomial regression along the time axis within each segment.
! kcode: 0 uniform (Savitzky-Golay / moving average), 1 Gaussian (band=sigma),
!        2 tricube (band=span in (0,1)).
! degree: 0 = local weighted mean, 1 = linear (loess), 2/3 = Savitzky-Golay.
! ---------------------------------------------------------------------------
subroutine smooth_locreg(x, nobs, nfeat, ntime, kcode, degree, band, nbrk, brk, y)
  implicit none
  integer, intent(in)  :: nobs, nfeat, ntime, kcode, degree, nbrk
  double precision, intent(in)  :: x(nobs, nfeat, ntime)
  double precision, intent(in)  :: band
  integer, intent(in)  :: brk(*)
  double precision, intent(out) :: y(nobs, nfeat, ntime)
  integer :: nseg, s, i, j, t, k, d, h, m, r, c
  integer, allocatable :: lo(:), hi(:)
  double precision :: aa(degree + 1, degree + 1), bb(degree + 1), solv(degree + 1)
  double precision :: w, u, spanh
  allocate(lo(nbrk + 1), hi(nbrk + 1))
  call build_segments(ntime, nbrk, brk, nseg, lo, hi)
  m = degree + 1
  do j = 1, nfeat
    do i = 1, nobs
      do s = 1, nseg
        if (hi(s) < lo(s)) cycle
        if (kcode == 2) then
          ! tricube span scales with the number of intervals in the segment
          spanh = band * dble(hi(s) - lo(s)) / 2.0d0
          h = max(1, int(floor(spanh)))
        else if (kcode == 1) then
          h = max(1, nint(3.0d0 * band))
        else
          h = max(0, nint(band))
        end if
        do t = lo(s), hi(s)
          aa = 0.0d0
          bb = 0.0d0
          do k = max(lo(s), t - h), min(hi(s), t + h)
            d = k - t
            if (kcode == 1) then
              u = dble(d) / band
              w = exp(-0.5d0 * u * u)
            else if (kcode == 2) then
              u = dble(d) / max(1.0d0, spanh)
              if (abs(u) <= 1.0d0) then
                w = (1.0d0 - abs(u) ** 3) ** 3
              else
                w = 0.0d0
              end if
            else
              w = 1.0d0
            end if
            do r = 1, m
              do c = 1, m
                aa(r, c) = aa(r, c) + w * dble(d) ** (r + c - 2)
              end do
              bb(r) = bb(r) + w * x(i, j, k) * dble(d) ** (r - 1)
            end do
          end do
          call solve_poly(aa, bb, m, solv)
          y(i, j, t) = solv(1)
        end do
      end do
    end do
  end do
  deallocate(lo, hi)
end subroutine smooth_locreg

! ---------------------------------------------------------------------------
! Exponentially weighted moving average. The recursion resets at each segment
! start so an intervention effect is preserved.
! ---------------------------------------------------------------------------
subroutine smooth_ewma(x, nobs, nfeat, ntime, alpha, nbrk, brk, y)
  implicit none
  integer, intent(in)  :: nobs, nfeat, ntime, nbrk
  double precision, intent(in)  :: x(nobs, nfeat, ntime)
  double precision, intent(in)  :: alpha
  integer, intent(in)  :: brk(*)
  double precision, intent(out) :: y(nobs, nfeat, ntime)
  integer :: nseg, s, i, j, t
  integer, allocatable :: lo(:), hi(:)
  double precision :: acc
  allocate(lo(nbrk + 1), hi(nbrk + 1))
  call build_segments(ntime, nbrk, brk, nseg, lo, hi)
  do j = 1, nfeat
    do i = 1, nobs
      do s = 1, nseg
        if (hi(s) < lo(s)) cycle
        acc = x(i, j, lo(s))
        y(i, j, lo(s)) = acc
        do t = lo(s) + 1, hi(s)
          acc = alpha * x(i, j, t) + (1.0d0 - alpha) * acc
          y(i, j, t) = acc
        end do
      end do
    end do
  end do
  deallocate(lo, hi)
end subroutine smooth_ewma

! ---------------------------------------------------------------------------
! Penalized cubic-spline-type smoother (second-difference penalty, Whittaker /
! Hodrick-Prescott style) per segment. Solves (I + lambda * D2' D2) z = y by
! banded (half-bandwidth 2) Cholesky, then y = z. lambda -> 0 interpolates,
! lambda -> inf forces a linear trend.
! ---------------------------------------------------------------------------
subroutine smooth_spline(x, nobs, nfeat, ntime, lambda, nbrk, brk, y)
  implicit none
  integer, intent(in)  :: nobs, nfeat, ntime, nbrk
  double precision, intent(in)  :: x(nobs, nfeat, ntime)
  double precision, intent(in)  :: lambda
  integer, intent(in)  :: brk(*)
  double precision, intent(out) :: y(nobs, nfeat, ntime)
  integer :: nseg, s, i, j, t, n
  integer, allocatable :: lo(:), hi(:)
  double precision, allocatable :: d0(:), d1(:), d2(:)
  double precision, allocatable :: ld0(:), ld1(:), ld2(:), rhs(:), ww(:), zz(:)
  double precision :: sm
  allocate(lo(nbrk + 1), hi(nbrk + 1))
  call build_segments(ntime, nbrk, brk, nseg, lo, hi)
  do j = 1, nfeat
    do i = 1, nobs
      do s = 1, nseg
        n = hi(s) - lo(s) + 1
        if (n < 2) then
          if (n == 1) y(i, j, lo(s)) = x(i, j, lo(s))
          cycle
        end if
        if (n == 2) then
          ! no second differences exist for a length-2 segment: identity
          y(i, j, lo(s)) = x(i, j, lo(s))
          y(i, j, lo(s) + 1) = x(i, j, lo(s) + 1)
          cycle
        end if
        allocate(d0(n), d1(n), d2(n), ld0(n), ld1(n), ld2(n), rhs(n), ww(n), zz(n))
        do t = 1, n
          rhs(t) = x(i, j, lo(s) + t - 1)
          d0(t) = 1.0d0 + lambda * 6.0d0
          d1(t) = -lambda * 4.0d0
          d2(t) = lambda
        end do
        d0(1) = 1.0d0 + lambda * 1.0d0
        d0(n) = 1.0d0 + lambda * 1.0d0
        d0(2) = 1.0d0 + lambda * 4.0d0
        if (n >= 4) d0(2) = 1.0d0 + lambda * 5.0d0
        if (n >= 4) d0(n - 1) = 1.0d0 + lambda * 5.0d0
        d1(1) = -lambda * 2.0d0
        d1(n - 1) = -lambda * 2.0d0
        ! banded Cholesky M = L L' (half-bandwidth 2); store ld0=diag, ld1=sub-1,
        ! ld2=sub-2. For row t: L[t,t-2]=d2(t-2)/L[t-2,t-2];
        ! L[t,t-1]=(d1(t-1)-L[t,t-2]*L[t-1,t-2])/L[t-1,t-1]; L[t,t]=sqrt(...).
        ld0 = 0.0d0; ld1 = 0.0d0; ld2 = 0.0d0
        do t = 1, n
          if (t >= 3) ld2(t) = d2(t - 2) / ld0(t - 2)
          if (t >= 2) then
            sm = d1(t - 1)
            if (t >= 3) sm = sm - ld2(t) * ld1(t - 1)
            ld1(t) = sm / ld0(t - 1)
          end if
          sm = d0(t)
          if (t >= 2) sm = sm - ld1(t) ** 2
          if (t >= 3) sm = sm - ld2(t) ** 2
          ld0(t) = sqrt(sm)
        end do
        ! forward solve L w = rhs
        ww(1) = rhs(1) / ld0(1)
        do t = 2, n
          sm = rhs(t) - ld1(t) * ww(t - 1)
          if (t >= 3) sm = sm - ld2(t) * ww(t - 2)
          ww(t) = sm / ld0(t)
        end do
        ! backward solve L' z = w
        zz(n) = ww(n) / ld0(n)
        do t = n - 1, 1, -1
          sm = ww(t) - ld1(t + 1) * zz(t + 1)
          if (t <= n - 2) sm = sm - ld2(t + 2) * zz(t + 2)
          zz(t) = sm / ld0(t)
        end do
        do t = 1, n
          y(i, j, lo(s) + t - 1) = zz(t)
        end do
        deallocate(d0, d1, d2, ld0, ld1, ld2, rhs, ww, zz)
      end do
    end do
  end do
  deallocate(lo, hi)
end subroutine smooth_spline

! ---------------------------------------------------------------------------
! Paper-depth Kalman family: hierarchical negative-binomial EKF / RTS.
! The latent state z is the log relative abundance (CLR-equivalent up to a
! per-sample constant), the expected count is mean = lib * exp(z) (library size
! times the softmax of the state), and the measurement variance follows the
! delta method, var = mean + mean^2 / phi, with per-feature dispersion phi.
! A hierarchical EM loop (shared per-feature innovation variance q and
! dispersion phi across all observations, ~80 iterations in the paper) is run
! first when n_em > 0; the process noise is optionally scaled by the sampling
! interval dt for irregular schedules. The output y is the smoothed expected
! count (lib * exp(z_smoothed)), keeping the count scale of the input.
! ---------------------------------------------------------------------------
subroutine smooth_kalman_nb(x, nobs, nfeat, ntime, lib, phi, q, mode, &
                            n_em, dt, use_dt, y)
  implicit none
  integer, intent(in)  :: nobs, nfeat, ntime, mode, n_em, use_dt
  double precision, intent(in)  :: x(nobs, nfeat, ntime)
  double precision, intent(in)  :: lib(nobs, ntime)
  double precision, intent(inout) :: phi(nfeat), q(nfeat)
  double precision, intent(in)  :: dt(ntime)
  double precision, intent(out) :: y(nobs, nfeat, ntime)
  double precision, parameter :: floor_ = -16.0d0, half = 0.5d0
  double precision, parameter :: minphi = 1.0d-4, minq = 1.0d-6, emtol = 1.0d-6
  double precision :: zsm(nobs, nfeat, ntime)
  double precision :: zf(ntime), pf(ntime), zp(ntime), pp(ntime), zs(ntime)
  double precision :: mu, h, s, k, gain, diffq, diffphi, acc
  double precision :: qnew, phinew, sumnum, sumden
  integer :: it, i, f, t
  ! Initial latent state = log relative abundance (pseudo-counted at zero).
  do f = 1, nfeat
    do i = 1, nobs
      do t = 1, ntime
        zsm(i, f, t) = log((x(i, f, t) + half) / lib(i, t))
      end do
    end do
  end do
  if (n_em > 0) then
    do it = 1, n_em
      ! E-step: RTS smoother for every (observation, feature) with current (q, phi).
      do f = 1, nfeat
        do i = 1, nobs
          zf(1) = zsm(i, f, 1)
          pf(1) = max(q(f), minq)
          zp(1) = zf(1); pp(1) = pf(1)
          do t = 2, ntime
            zp(t) = zf(t - 1)
            pp(t) = pf(t - 1) + q(f) * dble(merge(dt(t), 1.0d0, use_dt == 1))
            mu = lib(i, t) * exp(max(zp(t), floor_))
            h = mu
            s = h * h * pp(t) + mu + mu * mu / max(phi(f), minphi)
            k = pp(t) * h / s
            zf(t) = zp(t) + k * (x(i, f, t) - mu)
            pf(t) = (1.0d0 - k * h) * pp(t)
          end do
          zs(ntime) = zf(ntime)
          do t = ntime - 1, 1, -1
            gain = pf(t) / max(pp(t + 1), 1.0d-30)
            zs(t) = zf(t) + gain * (zs(t + 1) - zp(t + 1))
          end do
          do t = 1, ntime
            zsm(i, f, t) = max(zs(t), floor_)
          end do
        end do
      end do
      ! M-step: shared per-feature innovation variance and dispersion by moments.
      diffq = 0.0d0; diffphi = 0.0d0
      do f = 1, nfeat
        qnew = 0.0d0
        do i = 1, nobs
          do t = 2, ntime
            qnew = qnew + (zsm(i, f, t) - zsm(i, f, t - 1)) ** 2
          end do
        end do
        qnew = qnew / max(dble(nobs * max(ntime - 1, 1)), 1.0d0)
        qnew = max(qnew, minq)
        diffq = max(diffq, abs(qnew - q(f)) / max(q(f), 1.0d-30))
        q(f) = qnew
        sumnum = 0.0d0; sumden = 0.0d0
        do i = 1, nobs
          do t = 1, ntime
            mu = lib(i, t) * exp(zsm(i, f, t))
            if (mu <= 0.0d0) cycle
            acc = (x(i, f, t) - mu) ** 2 - mu
            if (acc > 0.0d0) then
              sumden = sumden + mu * mu
              sumnum = sumnum + mu * mu * mu / acc
            end if
          end do
        end do
        if (sumden > 0.0d0) then
          phinew = sumnum / sumden
          phinew = min(max(phinew, minphi), 1.0d3)
          diffphi = max(diffphi, abs(phinew - phi(f)) / max(phi(f), 1.0d-30))
          phi(f) = phinew
        end if
      end do
      if (max(diffq, diffphi) < emtol) exit
    end do
  end if
  ! Final pass: EKF forward (mode 0) or RTS (mode 1) with the estimated (q, phi).
  do f = 1, nfeat
    do i = 1, nobs
      zf(1) = zsm(i, f, 1)
      pf(1) = max(q(f), minq)
      zp(1) = zf(1); pp(1) = pf(1)
      do t = 2, ntime
        zp(t) = zf(t - 1)
        pp(t) = pf(t - 1) + q(f) * dble(merge(dt(t), 1.0d0, use_dt == 1))
        mu = lib(i, t) * exp(max(zp(t), floor_))
        h = mu
        s = h * h * pp(t) + mu + mu * mu / max(phi(f), minphi)
        k = pp(t) * h / s
        zf(t) = zp(t) + k * (x(i, f, t) - mu)
        pf(t) = (1.0d0 - k * h) * pp(t)
      end do
      if (mode == 0) then
        do t = 1, ntime
          y(i, f, t) = lib(i, t) * exp(max(zf(t), floor_))
        end do
      else
        zs(ntime) = zf(ntime)
        do t = ntime - 1, 1, -1
          gain = pf(t) / max(pp(t + 1), 1.0d-30)
          zs(t) = zf(t) + gain * (zs(t + 1) - zp(t + 1))
        end do
        do t = 1, ntime
          y(i, f, t) = lib(i, t) * exp(max(zs(t), floor_))
        end do
      end if
    end do
  end do
end subroutine smooth_kalman_nb

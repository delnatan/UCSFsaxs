subroutine autorg(q,y,Nmin,qRgmax,Npts,I0,opt_Rg,opt_qRg,gi,gf)
    integer Npts,Nmin,blocksize,ii,ei,opt_block,opt_ii,Nthird
    real(8), dimension(Npts) :: q,y,qsq,logy
    real(8) :: qRgmax,m,b,r,merr,berr,min_merr,rg,qrg
    real(8), intent(out) :: opt_Rg,I0,opt_qRg
    integer, intent(out) :: gi,gf
    ! output slope, y-intercept, optimal q*Rg, guinier_begin, guinier_end
    c = 1
    min_merr  = 100.0
    blocksize = Nmin
    opt_block = Nmin
    ! explore one-third of dataset
    Nthird = Npts/3
    ii = 1

    qsq = q**2
    logy= log(y)

    do while (blocksize .lt. Nthird)
        ii = 1
!        write (*,102) blocksize
!102     format("Working with block-size ", i3,".")
        do while (ii.lt.(Nthird-blocksize))
            ei = ii + blocksize
            call linregress(qsq(ii:ei), logy(ii:ei), blocksize, m, b, r, merr, berr)
            rg  = sqrt(-3.0 * m)
            qrg = q(ei) * rg
            if (merr .lt. min_merr .and. qrg.le.qRgmax .and. m.lt.0) then
                min_merr = merr
                opt_block = blocksize
                opt_ii    = ii
                opt_qRg   = qrg
                opt_Rg    = rg
                I0        = exp(b)
            end if
            ii  = ii + 1
        end do
        blocksize = blocksize + 1
    end do

    gi = opt_ii
    gf = opt_ii + opt_block

    write(*,101) opt_ii, opt_ii+opt_block, opt_qRg, opt_Rg
101 format('Rg (',i3,' - ',i3,')',1x,' qRg=',f6.3,', Rg=',f8.3)

end subroutine autorg


subroutine linregress(x, y, N, slope, intercept, r, stderr_slope, stderr_intercept)
! Calculates simple linear regression of (x, y)
    integer :: N
    real(8), dimension(N), intent(in)  :: x,y
    real(8), intent(out) :: slope 
    real(8), intent(out) :: intercept 
    real(8), intent(out) :: r 
    real(8), intent(out) :: stderr_slope
    real(8), intent(out) :: stderr_intercept
    real(8) :: xmean, ymean, varx, covxy, vary, r_den, mse


    N = size(x)
    xmean = sum(x)/N
    ymean = sum(y)/N
    varx = dot_product(x-xmean, x-xmean)
    covxy = dot_product(x-xmean, y-ymean)
    vary = dot_product(y-ymean, y-ymean)

    slope = covxy / varx
    intercept = ymean - slope*xmean

    r_den = sqrt(varx * vary)
    if (abs(r_den) < tiny(1.0)) then
        r = 0
    else
        r = covxy / r_den
        ! Normalize to [-1, 1] in case of numerical error propagation
        if (r > 1) then
            r = 1
        else if (r < -1) then
            r = -1
        end if
    endif
    ! 'mse' is a mean square error (the sum of squared residuals divided by number
    ! of model parameters), which can be calculated directly as:
    !   mse = sum((y-(slope*x+intercept))**2) / (N-2)
    ! But equivalently it can also be calculated in a faster way as:
    mse = (1-r**2) * vary / (N-2)
    stderr_slope = sqrt(mse / varx)
    stderr_intercept = sqrt(mse * (1.0/N + xmean**2/varx))

end subroutine linregress

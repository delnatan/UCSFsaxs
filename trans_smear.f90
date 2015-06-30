
subroutine transc(q,r,y,Wy,Nq,Nr,Ny,A)
    ! smearing correction taken into account 
    ! beam parameters are y and Wy
    real(8), parameter :: pi = acos(-1.0)

    integer :: Nq,Nr,Ny,i,j,k
    real(8), dimension(Nq), intent(in) :: q
    real(8), dimension(Nr), intent(in) :: r
    real(8), dimension(Ny), intent(in) :: y, Wy
    real(8), dimension(Nq,Nr), intent(out) :: A
    real(8) :: qy, aij, dr, dy

    dr = r(2)-r(1)
    dy = y(2)-y(1)

    do i=1,Nq
        do j=1,Nr
            aij = 0.0
            do k=1,Ny
                qy = sqrt(q(i)**2 +y(k)**2) * r(j)
                if (qy<1e-15) then
                    aij = aij + Wy(k) * 1.0 * dy
                else
                    aij = aij + Wy(k) * sin(qy)/qy * dy
                end if
            end do
            A(i,j) = 4.0*pi*2.0*aij*dr
        end do
    end do

end subroutine transc

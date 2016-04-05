!%fortran
!F90
    subroutine omegaf_vec(Cpp,xg,yg,ng,x_meas,y_meas,n_meas,L,dist)
        integer ng,n_meas,i,j
        real :: L
        real :: xg(ng),yg(ng),x_meas(n_meas),y_meas(n_meas),Cpp(ng,n_meas)
        real :: dist(ng,n_meas)
        real a
        real b(ng,n_meas),boa(ng,n_meas)
        real c(ng,n_meas)
!f2py   intent(in)     :: L
!f2py   intent(in)     :: xg,yg,x_meas,y_meas
!f2py   intent(in,out) :: Cpp
!f2py   intent(out)    :: dist

        if(L.eq.0) then
         return
        end if

        ! for very large memory cases, the following script might be better,
        ! % as it avoids using repmat
        do i=1,ng
          do j=1,n_meas
              dist(i,j)=sqrt( (xg(i)-x_meas(j))**2+(yg(i)-y_meas(j))**2 )
          end do
        end do

        a = sqrt(10.0/3.0)*L
        b = dist
        boa=b/a
        c = 0*Cpp
        do i=1,ng
          do j=1,n_meas
               if ((0 .le. b(i,j)) .and. (b(i,j) .le. a)) then
                   c(i,j) =                 &
                         -(1/4)*boa(i,j)**5 &
                         +(1/2)*boa(i,j)**4 &
                         +(5/8)*boa(i,j)**3 &
                         -(5/3)*boa(i,j)**2 &
                         +1
                end if
          end do
        end do

        do i=1,ng
          do j=1,n_meas
               if ((a .le. b(i,j)) .and. (b(i,j) .le. 2 * a)) then
                   c(i,j) =               &
                       -(1/4)*boa(i,j)**5 &
                       +(1/2)*boa(i,j)**4 &
                       +(5/8)*boa(i,j)**3 &
                       -(5/3)*boa(i,j)**2 &
                       +1
                end if

          end do
        end do

        Cpp=Cpp*c

    continue
    end

!%fortran
!F90
    subroutine distg(xg,yg,ng,x_meas,y_meas,n_meas,dist)
        integer ng,n_meas,i,j
        real :: xg(ng),yg(ng),x_meas(n_meas),y_meas(n_meas)
        real :: dist(ng,n_meas)
!f2py   intent(in)     :: xg,yg,x_meas,y_meas
!f2py   intent(out)    :: dist

        do i=1,ng
          do j=1,n_meas
              dist(i,j)=sqrt( (xg(i)-x_meas(j))**2+(yg(i)-y_meas(j))**2 )
          end do
        end do
    continue
    end


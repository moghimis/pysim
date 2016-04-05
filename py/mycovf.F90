!%fortran
!F90

    subroutine mycovf(a,na,b,nb,nn,cab)
        integer  ::  na,nb,nn,i,j
        real     ::  a(na,nn), b(nb,nn)
        real     ::  da(na,nn),db(nb,nn)
        real     ::  am(na,nn), bm(nb,nn)
        real     ::  cab(na,nb)
!f2py   intent(in)     :: a,b,na,nb,nn
!f2py   intent(out) :: cab

        !function Cab = myCov(a,b)
        !
        ! Cab = myCov(a,b)
        !
        ! Computes sample covariance 'a' and 'b'.  Inputs should have dimensions
        ! as follows...
        !
        !    size(a) = [ M, N ]
        !    size(b) = [ P, N ]
        !
        ! That is, each matrix consists of N samples of an Mx1 (or Px1) vector.
        ! Note M and P need not be equal (unlike the builtin matlab function cov.m).
        ! Output 'Cab' is the sample covariance, which has dimensions MxP.
        !if(size(a,2) ~= size(b,2))
        !  error('inputs must have same number of samples (2nd dimension)')
        !end
        !
        !da=(a-repmat(nanmean(a,2),1,size(a,2)));
        !db=(b-repmat(nanmean(b,2),1,size(b,2)));
        !N=size(a,2);
        !Cab = (1/(N-1))*da*db';

        do i=1,na
           am(i,j) = sum(a(i,:))/size(a(i,:))
        end do

        do i=1,nb
           bm(i,j) = sum(b(i,:))/size(b(i,:))
        end do

        do i=1,na
            do j=1,nn
                da(i,j) = a(i,j) - am(i,j)
            end do
        end do

        do i=1,nb
            do j=1,nn
                db(i,j) = b(i,j) - bm(i,j)
            end do
        end do

        da=1
        db=1
        cab=(1/(nn-1))* matmul(da,transpose(db))
    continue
    end



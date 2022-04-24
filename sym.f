        module sym
        implicit none
        save
        real*8, intent(in) :: Cb(:,:)
        contains


        subroutine rotation(Cb,natom)
        implicit none
        integer :: natom
        allocate(Cb(natom,3))
        do n=6,1
        pi=4.d0*DATAN(1.d0)
        Cn(1,1)=dcos(2*pi/n)
        Cn(1,2)=dsin(2*pi/n)
        Cn(2,2)=dcos(2*pi/n)
        Cn(2,1)=-dsin(2*pi/n)
        Cn(3,3)=1.d0
        R(n)=Cn(i,j)
        enddo




        end subroutine

        end module


        program uno
        implicit none
        integer :: i, j, natom,n,k
        real*8, allocatable :: A(:),M(:)
        real*8, allocatable :: C(:,:),Ia(:,:),II(:,:)
        real*8, allocatable :: Ca(:,:,:),Pa(:,:,:)
        real*8, allocatable :: Cn(:,:,:)!JUSTincase
        real*8 :: pi,acum
        character(LEN=2) :: label1,label2
        character(LEN=2),allocatable :: element(:) 
        pi=4.d0*DATAN(1.d0)
        write(6,*)pi
        open(1,file='input.mol')
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)natom
        allocate(A(3),M(3),C(natom,3),element(natom))
        allocate(Pa(2,natom,3),Ia(natom,3),Ca(6,natom,3))
        allocate(Cn(6,3,3)) !JUSTincase
        if(natom.eq.2)then
        read(1,*)A(1),A(2),A(3),label1
        read(1,*)M(1),M(2),M(3),label2
                if(label1.eq.label2)then
                        write(6,*)'The molecule belongs Dinfh'
                else
                       write(6,*)'The molecule belongs Cinfv'
                endif
        endif
        do i=1,natom
        read(1,*)(C(i,j),j=1,3),element(i)
        enddo
        write(6,*)'ok',natom,C(3,1)
        !aquí tenemos que ordenar la matriz
        call inversion(C,natom,Ia)
        write(6,*)C(3,1),Ia(3,1),Ia(1,1) 
        !def rotation int------------------------------
        do n=1,6
        Cn(n,i,j)=0.d0
        Cn(n,1,1)=dcos(2*pi/n)
        Cn(n,1,2)=-dsin(2*pi/n)
        Cn(n,2,2)=dcos(2*pi/n)
        Cn(n,2,1)=dsin(2*pi/n)
        Cn(n,3,3)=1.d0
        enddo
        !matrix multiplication-------------------------
        do n=1,6
          do i=1,natom
          do j=1,3
          acum=0.d0
            do k=1,3
            acum=acum+C(i,k)*Cn(n,k,j)
            enddo
            Ca(n,i,j)=acum
          enddo
          enddo
         enddo
        write(6,*)C(3,1),Ca(1,3,1),C(3,1),Ca(3,3,1),C(3,3),Ca(3,3,3)
        end program

        subroutine rotation(C,natom,Ca)
        implicit none
        real*8 :: pi, acum
        real*8, intent(in) :: C(natom,3)
        real*8, intent(out) :: Ca(6,natom,3)
        real*8, allocatable :: Cn(:,:,:)
        integer :: natom,i,j,n,k
        allocate(Cn(6,3,3))
        pi=4.d0*DATAN(1.d0)
        do n=1,6
        Cn(n,i,j)=0.d0
        Cn(n,1,1)=dcos(2*pi/n)
        Cn(n,1,2)=-dsin(2*pi/n)
        Cn(n,2,2)=dcos(2*pi/n)
        Cn(n,2,1)=dsin(2*pi/n)
        Cn(n,3,3)=1.d0
        enddo
        do n=1,6
          do i=1,natom
          do j=1,3
          acum=0.d0
            do k=1,3
            acum=acum+C(i,k)*Cn(n,k,j)
            enddo
            Ca(n,i,j)=acum
          enddo
          enddo
        enddo
        end subroutine 

        subroutine inversion2(C,natom,Ia)
        implicit none
        real*8, intent(in) :: C(natom,3)
        real*8, intent(out) :: Ia(natom,3)
        integer :: natom,i,j,n,k
        !II(i,j)=0.d0
          do i=1,natom
          do j=1,3
          Ia(i,j)=-C(i,j)
          enddo
          enddo
        end subroutine

        subroutine inversion(C,natom,Ia)
        implicit none
        real*8 :: acum
        real*8, intent(in) :: C(natom,3)
        real*8, intent(out) :: Ia(natom,3)
        real*8, allocatable :: II(:,:)
        integer :: natom,i,j,n,k
        allocate(II(3,3))
         do i=1,natom
         do j=1,3
          if(i.eq.j)then
          II(i,j)=-1.d0
          else
          II(i,j)=0.d0
          endif
          !acum=0.d0
          ! do k=1,3
           ! acum=acum+C(i,k)*II(k,j)
           !enddo
          enddo
          enddo
           II(i,j)=Ia(i,j)
           return
        end subroutine 

        subroutine planes(C,natom,Pa)
        implicit none
        real*8 :: acum
        real*8, intent(in) :: C(natom,3)
        real*8, intent(out) :: Pa(3,natom,3)
        real*8, allocatable :: Pn(:,:,:)
        integer :: natom,i,j,n,k
        allocate(Pn(3,3,3))
        !Pn(1)=sigma_yz Pn(2)=sigma_xz Pn(3)=sigma_zy
        Pn(n,i,j)=0.d0
        Pn(1,1,1)=-1
        Pn(1,2,2)=1
        Pn(1,3,3)=1
        Pn(2,1,1)=1
        Pn(2,2,2)=-1
        Pn(2,3,3)=1
        Pn(3,1,1)=1
        Pn(3,2,2)=1
        Pn(3,3,3)=-1
        do n=1,3
          do i=1,natom
          do j=1,3
          acum=0.d0
            do k=1,3
            acum=acum+C(i,k)*Pn(n,k,j)
            enddo
          Pa(n,i,j)=acum
          enddo
          enddo
        enddo
        end subroutine

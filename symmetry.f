        program symmetry
        implicit none
        integer :: i,j,natom,n,k,l,p,q,s,r,v,x
        integer, allocatable :: axis(:),order(:)
        real*8, allocatable :: A(:),M(:)
        real*8, allocatable :: C(:,:),Ia(:,:)
        real*8, allocatable :: Ca(:,:,:,:),Pa(:,:,:)
        real*8 :: pi,indx,indy,indz,subx,suby,subz,diffx,diffy
        character(LEN=2) :: label1,label2,label
        character(LEN=2),allocatable :: element(:) 
        pi=4.d0*DATAN(1.d0)
        !open(1,file='inputBEN.xyz') !benzene
        !open(1,file='inputH2O.xyz') !water
        !open(1,file='inputCYCLOB.xyz') !cyclobutane
        !open(1,file='inputDICLORO.xyz') !C1
        open(1,file='inputTOLUENO.mol') !Cs
        read(1,*)natom
        read(1,*)
        allocate(A(3),M(3),C(natom,3),element(natom))
        allocate(Pa(2,natom,3),Ia(natom,3),Ca(3,6,natom,3))
        !---Diatomic molecules ONLY -------------------------------------
        if(natom.eq.2)then
        read(1,*)A(1),A(2),A(3),label1
        read(1,*)M(1),M(2),M(3),label2
                if(label1.eq.label2)then
                        write(6,*)'The molecule belongs Dinfh'
                else
                       write(6,*)'The molecule belongs Cinfv'
                endif
        endif
        !-----Coordinates Reading----------------------------------------
        do i=1,natom
        read(1,*)element(i),(C(i,j),j=1,3)
        write(6,*)(C(i,j),j=1,3)
        enddo
        
        call rotation(C,natom,Ca)
        !MATRIX Ca(i,j,k,l) where i=rotation respect the axe 1:x 2:y 3:z
                                 !j=order of the rotation [1-6]
                                 !k=rows of atoms 
                                 !l=position x,y and z. 

        !------ Sorting ------------------------------------------------
        indx=2000
        indy=2000
        indz=2000
        do i=1,3
        do j=1,6
        do q=1,natom
           do p=1,natom
           
        if(p.ge.q)then
              diffx=indx-Ca(i,j,p,1)
              diffy=indy-Ca(i,j,p,2)
              if(diffx .gt. 0.0001)then
                    indx=Ca(i,j,p,1)
                    indy=Ca(i,j,p,2)
                    indz=Ca(i,j,p,3)
                    label=element(p)
                    s=p
              endif

              if(diffx .lt. -0.0001)then
                      continue
              endif

              if(diffx .lt. 0.0001 .and. diffx .gt. -0.0001             &
     &              .and.Ca(i,j,p,2).ne.indy)then
              indx=Ca(i,j,p,1)
               if(diffy .gt. 0.0001)then
                      indy=Ca(i,j,p,2)
                      indz=Ca(i,j,p,3)
                      label=element(p)
                      s=p
               endif
              endif
              if(diffx .lt. 0.0001 .and. diffx .gt. -0.0001             &
     &           .and. diffy .lt. 0.0001 .and. diffy .gt. -0.0001)then
                      indx=Ca(i,j,p,1)
                      indy=Ca(i,j,p,2)
                      if(Ca(i,j,p,3).lt.indz)then
                              indz=Ca(i,j,p,3)
                              label=element(p)
                              s=p
                      endif
              endif         
        endif
        enddo
        Ca(i,j,s,1)=Ca(i,j,q,1)
        Ca(i,j,s,2)=Ca(i,j,q,2)
        Ca(i,j,s,3)=Ca(i,j,q,3)
        element(s)=element(q)
        Ca(i,j,q,1)=indx
        Ca(i,j,q,2)=indy
        Ca(i,j,q,3)=indz
        element(q)=label
        indx=2000
        indy=2000
        indz=2000
        enddo
        enddo
        enddo
        !write(6,*)'Original matrix sorted:'
        do i=1,natom
        !write(6,*)(Ca(3,1,i,j),j=1,3),element(i)
        enddo
        !write(6,*)'Rotated by C^2(z) and sorted:'
        do i=1,natom
        !write(6,*)(Ca(3,2,i,j),j=1,3),element(i)
        enddo
        !write(6,*)'Rotated by C^2(x) and sorted:'
        do i=1,natom
        !write(6,*)(Ca(1,2,i,j),j=1,3),element(i)
        enddo
        !write(6,*)'Rotated by C^2(z) and sorted:'
        do i=1,natom
        !write(6,*)(Ca(3,2,i,j),j=1,3),element(i)
        enddo
        call compare(Ca,natom)
        open(88,file='count.dat')
        open(89,file='test.dat')
        n=0
        do while(.true.)
        read(88,*,end=10)
        n=n+1
        enddo
10      continue
        write(6,*)n
!Low symmetry groups----------------------------
        if(n.eq.0)then
                call inversion(C,natom,Ia)
        do i=1,natom
        write(6,*)
        enddo
        r=0
        do i=1,natom
        do j=1,natom
        if(C(i,1).eq.Ia(j,1).and.C(i,2).eq.Ia(j,2).and.C(i,3).eq.       &
     & Ia(j,3))then
                r=r+1
        endif
        enddo
        enddo 
        if(r.gt.2)then
        write(6,*)'this molecule belongs to the C_i symmetry group'
        else

        write(6,*)'this molecule belongs to the C_1 symmetry group'    
        endif

        endif
        
!--------------------------------------------------------
        call planes(C,natom,Pa)
        write(6,*)'original'       
        do i=1,natom
        write(6,*)(C(i,j),j=1,3)       
        enddo
        do n=1,3
        write(6,*)n
        do i=1,natom
        write(6,*)(Pa(n,i,j),j=1,3)
        enddo
        write(6,*)
        enddo
        
        allocate(axis(n),order(n))
        do i=1,n
                read(89,*)axis(i),order(i)
        enddo
        x=maxval(axis)
        do i=1,n
            if(axis(i).eq.x)then
                  !  write(6,*)axis(i),order(i)
            endif
        enddo


               
        
!-----Calling -------------------------------------------------
        call inversion(C,natom,Ia)
        do i=1,natom
         write(6,*)(C(i,j),j=1,3)
         write(6,*)(Ia(i,j),j=1,3)
         write(6,*)
        enddo


        
        end program

        !**************************************************************!



        subroutine rotation(C,natom,Ca)
        implicit none
        real*8 :: pi, acum
        real*8, intent(in) :: C(natom,3)
        real*8, intent(out) :: Ca(3,6,natom,3)
        real*8, allocatable :: Cn(:,:,:,:)
        integer :: i,j,n,k,natom,l
        allocate(Cn(3,6,3,3))
        pi=4.d0*DATAN(1.d0)
        Cn=0.d0
        !def rotacion x-----------------
        do n=1,6
        Cn(1,n,1,1)=1.d0
        Cn(1,n,2,3)=dsin(2*pi/n)
        Cn(1,n,2,2)=dcos(2*pi/n)
        Cn(1,n,3,2)=-dsin(2*pi/n)
        Cn(1,n,3,3)=dcos(2*pi/n)
        enddo
        !defrotation y-----------------
        do n=1,6
        Cn(2,n,1,1)=dcos(2*pi/n)
        Cn(2,n,1,3)=-dsin(2*pi/n)
        Cn(2,n,2,2)=1.d0
        Cn(2,n,3,1)=dsin(2*pi/n)
        Cn(2,n,3,3)=dcos(2*pi/n)
        enddo
        !def rotacion z----------------------
        do n=1,6
        Cn(3,n,1,1)=dcos(2*pi/n)
        Cn(3,n,1,2)=-dsin(2*pi/n)
        Cn(3,n,2,2)=dcos(2*pi/n)
        Cn(3,n,2,1)=dsin(2*pi/n)
        Cn(3,n,3,3)=1.d0
        enddo
        !matrix multiplication
        do l=1,3
        do n=1,6
          do i=1,natom
          do j=1,3
          acum=0.d0
            do k=1,3
            acum=acum+C(i,k)*Cn(l,n,k,j)
            enddo
            Ca(l,n,i,j)=acum
          enddo
          enddo
         enddo
        enddo
        end subroutine 
        
        
        subroutine planes(C,natom,Pa)
        implicit none
        real*8 :: pi, acum
        real*8, intent(in) :: C(natom,3)
        real*8, intent(out) :: Pa(3,natom,3)
        real*8, allocatable :: Pn(:,:,:)
        integer :: i,j,n,k,natom
        allocate(Pn(3,3,3))
        Pn=0.d0
        !Pn(1,i,j)-yz|Pn(2,i,j)-xz|Pn(3,i,j)-xy| 
        do n=1,3
        do i=1,3
        do j=1,3
          if(i.eq.j)then
                  Pn(n,i,j)=1.d0
          endif
        Pn(1,1,1)=-1.d0
        Pn(2,2,2)=-1.d0
        Pn(3,3,3)=-1.d0
        enddo
        enddo
        enddo
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
        

        subroutine inversion(C,natom,Ia)
        implicit none
        real*8 :: acum
        real*8, intent(in) :: C(natom,3)
        real*8, intent(out) :: Ia(natom,3)
        real*8, allocatable :: II(:,:)
        integer :: natom,i,j,n,k
        allocate(II(3,3))
         II=0.d0
         do i=1,natom
         do j=1,3
           if(i.eq.j)then
             II(i,j)=-1.d0
           endif
           acum=0.d0
           do k=1,3
             acum=acum+C(i,k)*II(k,j)
           enddo
            Ia(i,j)=acum
         enddo
         enddo
        end subroutine 

        subroutine compare(Ca,natom)
        !----Compare C(i,1,k,l) with C(i,n,k,l)----------------------
        !Compara todas las filas y escribe el orden de la rotación(n) y
        !la dirección
        implicit none
        integer :: r,n,l,natom,i
        real*8,intent(in) :: Ca(3,6,natom,3)
        real*8 :: subx,suby,subz
        open(45,file='count.dat')
        open(46,file='test.dat')
        r=0
        do l=1,3     !direction
        do n=6,2,-1  !axis
        do i=1,natom !rows
         subx=Ca(l,1,i,1)-Ca(l,n,i,1)
         subx=abs(subx)
         suby=Ca(l,1,i,2)-Ca(l,n,i,2)
         suby=abs(suby)
         subz=Ca(l,1,i,3)-Ca(l,n,i,3)
         subz=abs(subz)
 
         if(subx .lt. 0.001 .and. suby .lt. 0.001 .and.                 &
     &      subz .lt. 0.001)then  
            r=r+1          
         
         endif
        enddo
        if(r .eq. natom)then
          write(45,*)n,l
          write(46,*)n,l
        end if
        
           r=0
        enddo
        enddo
        close(45)
        close(46)
        end subroutine 

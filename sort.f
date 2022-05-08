        program uno
        implicit none
        integer :: i, j, natom,p,q,s
        real*8, allocatable :: C(:,:)
        real*8 ::  indx,indy,indz
        character(LEN=2) :: label
        character(LEN=2),allocatable :: element(:) 
        open(1,file='inputBEN.xyz')
        read(1,*)natom
        read(1,*)
        allocate(C(natom,3),element(natom))
       
        write(6,*)'Sin ordenar'
        do i=1,natom
        read(1,*)element(i),(C(i,j),j=1,3)
        write(6,*)element(i),(C(i,j),j=1,3)
        enddo
        indx=2000
        indy=2000
        indz=2000
        write(6,*)'ordenada, i guess'
        do q=1,natom
                   do p=1,natom
                   
                   
                   IF (p .ge. q) THEN
                   
                      IF (C(p,1) .lt. indx) THEN
                          indx=C(p,1)
                          indy=C(p,2)
                          indz=C(p,3)
                          label=element(p)
                          s=p
                      END IF 
                      
                      
                      IF (C(p,1) .gt. indx) THEN
                          continue
                      END IF 
                      
                      IF (C(p,1) .eq. indx .and. C(p,2) .ne. indy) THEN
                          indx=C(p,1)
                          IF (C(p,2) .lt. indy) THEN
                              indy=C(p,2)
                              indz=C(p,3)
                              label=element(p)
                              s=p
                          END IF
                      
                      END IF
                      
                      IF (C(p,1) .eq. indx .and. C(p,2) .eq. indy) THEN
                          indx=C(p,1)
                          indy=C(p,2)
                          IF (C(p,3) .lt. indz) THEN
                              indz=C(p,3)
                              label=element(p)
                              s=p
                          END IF
                      
                      END IF
                  END IF
                  
                       
                     
                       
                  enddo
                  
               C(s,1)=C(q,1)
               C(s,2)=C(q,2)
               C(s,3)=C(q,3)
               element(s)=element(q)
               C(q,1)=indx
               C(q,2)=indy
               C(q,3)=indz
               element(q)=label
               indx=2000
               indy=2000
               indz=2000
               write(6,*)(C(q,j),j=1,3),element(q)
        enddo
        
        end program


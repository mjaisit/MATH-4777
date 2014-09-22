c**********************************************************************
c     matmat.f - matrix - matrix multiply, simple self-scheduling version
c************************************************************************
      program main
      include "mpif.h"

      integer MAX_AROWS, MAX_ACOLS, MAX_BCOLS 
      parameter (MAX_AROWS = 128, MAX_ACOLS = 128, MAX_BCOLS = 128)
      double precision a(MAX_AROWS,MAX_ACOLS), b(MAX_ACOLS,MAX_BCOLS)
      double precision c(MAX_AROWS,MAX_BCOLS)
      double precision buffer(MAX_ACOLS), ans(MAX_ACOLS)

      integer myid, master, numprocs, ierr, status(MPI_STATUS_SIZE)
      integer i, j, numsent, sender
      integer anstype, row, arows, acols, brows, bcols, crows, ccols

      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      print *, "Process ", myid, " of ", numprocs, " is alive"

      master = 0
      arows  = 128
      acols  = 128
      brows  = 128
      bcols  = 128
      crows  = arows
      ccols  = bcols

      if ( myid .eq. master ) then
c        master initializes and then dispatches
c        initialize a and b
         do 11 i = 1,acols
            do 10 j = 1,arows
               a(j,i) = i
 10         continue
 11      continue
         do 21 i = 1,bcols
            do 20 j = 1,brows
               b(j,i) = i
 20         continue
 21      continue

         numsent = 0
         
c        send b to each other process
         do 25 i = 1,bcols
         call MPI_BCAST(b(1,i), brows, MPI_DOUBLE_PRECISION, master, 
     &           MPI_COMM_WORLD, ierr)
 25      continue

c        send a row of a to each other process; tag with row number
         do 40 i = 1,numprocs-1
            do 30 j = 1,acols
               buffer(j) = a(i,j)
 30         continue
            call MPI_SEND(buffer, acols, MPI_DOUBLE_PRECISION, i, 
     &           i, MPI_COMM_WORLD, ierr)
            numsent = numsent+1
 40      continue
         
         do 70 i = 1,crows
            call MPI_RECV(ans, ccols, MPI_DOUBLE_PRECISION, 
     &                MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
     &            status, ierr)
            sender     = status(MPI_SOURCE)
            anstype    = status(MPI_TAG)
            do 45 j = 1,ccols
               c(anstype,j) = ans(j)
 45         continue

            if (numsent .lt. arows) then
               do 50 j = 1,acols
                  buffer(j) = a(numsent+1,j)
 50            continue
               call MPI_SEND(buffer, acols, MPI_DOUBLE_PRECISION, 
     &             sender, numsent+1, MPI_COMM_WORLD, ierr)
               numsent = numsent+1
            else
               call MPI_SEND(1.0, 1, MPI_DOUBLE_PRECISION, sender, 
     &             0, MPI_COMM_WORLD, ierr)
            endif
 70      continue
         
c        print out the answer
         do 80 i = 1,crows
            do 78 j = 1,ccols
               print *, "c(", i, ",", j, ") = ", c(i,j)
 78         continue
 80      continue

      else
c     slaves receive B, then compute rows of C until done message
      do 85 i = 1,bcols
         call MPI_BCAST(b(1,i), brows, MPI_DOUBLE_PRECISION, 
     &                   master, MPI_COMM_WORLD, ierr)
 85   continue
 90   call MPI_RECV(buffer, acols, MPI_DOUBLE_PRECISION, master, 
     &               MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      if (status(MPI_TAG) .eq. 0) then
         go to 200
      else
         row = status(MPI_TAG)
         do 100 i = 1,bcols
            ans(i) = 0.0
            do 95 j = 1,acols
               ans(i) = ans(i) + buffer(j)*b(j,i)
 95         continue
 100     continue
         call MPI_SEND(ans, bcols, MPI_DOUBLE_PRECISION, master, 
     &                  row, MPI_COMM_WORLD, ierr)
         go to 90
      endif
 200  continue
      endif

      call MPI_FINALIZE(ierr)
      stop
      end




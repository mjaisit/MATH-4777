c************************************************************************
c	Michael Anderson 
c	MATH 4777
c	HW 3
c	Martrix Matrix Multiplication
c	This program generates two random matrices A and B then computes
c	the product in parallel through the use of MPI
c************************************************************************
      program MM
      include "mpif.h"
	
c 	  First initialize the needed variables and set-up the matrices to the correct size
      integer A_rows, A_col, B_col 
      parameter (A_rows = 128, A_col = 128, B_col = 128)
      double precision a(A_rows,A_col), b(A_col,B_col),
		& c(A_rows,B_col)
		
c	  Initialize other needed variables for the operations	
	  integer anstype, row, arows, acols, brows, bcols, crows, ccols
      double precision temp(A_col), ans(A_col)
	  integer i, j
	  
c	  Initialize variables to use in open MPI
      integer master, numprocs, ierr, myid, status(MPI_STATUS_SIZE)
      integer numsent, sender
      

c	  Start MPI and throw error if not successful
      call MPI_INIT( ierr )
	  if (ierr .ne. MPI_SUCCESS) then
		print *,'Error starting MPI program. Terminating.'
		call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
	  end if
c	  Determine what the processor number is
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
c	  Determine how many processors there are and print each processors id
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      print *, "Process ", myid, " of ", numprocs, " is alive"

c	  Set the non parameter variables
      master = 0
      arows  = 128
      acols  = 128
      brows  = 128
      bcols  = 128
      crows  = arows
      ccols  = bcols


c	  Let the master set up the matrices A and B
c	  Each element in the matrix is going to set to the value of the column index
      if ( myid .eq. master ) then
         do 11 j = 1,acols
            do 10 i = 1,arows
               a(i,j) = j
 10         continue
 11      continue
         do 21 j = 1,bcols
            do 20 i = 1,brows
               b(i,j) = j
 20         continue
 21      continue


         numsent = 0


c        Send matrix B to the all other processors
         do 25 i = 1, bcols
         call MPI_BCAST(b(1,i), brows, MPI_DOUBLE_PRECISION, master, 
     &           MPI_COMM_WORLD, ierr)
 25      continue


c        send a row of a to each other process; tag with row number
         do 40 i = 1,numprocs-1
            do 30 j = 1,acols
               temp(j) = a(i,j)
 30         continue
            call MPI_SEND(temp, acols, MPI_DOUBLE_PRECISION, i, 
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
                  temp(j) = a(numsent+1,j)
 50            continue
               call MPI_SEND(temp, acols, MPI_DOUBLE_PRECISION, 
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
 90   call MPI_RECV(temp, acols, MPI_DOUBLE_PRECISION, master, 
     &               MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      if (status(MPI_TAG) .eq. 0) then
         go to 200
      else
         row = status(MPI_TAG)
         do 100 i = 1,bcols
            ans(i) = 0.0
            do 95 j = 1,acols
               ans(i) = ans(i) + temp(j)*b(j,i)
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








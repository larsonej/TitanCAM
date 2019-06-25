      program main
c
c
c  @(#) chekres.f  McKie  Oct-1995
c  This program checks the restart input and output source
c  code for the aerosol model to ensure that the same
c  variables are being input & output as are declared in
c  the global common blocks.
c
c  Modified Jun-1996:
c   Now uses template files for the restart input & output subroutines,
c   and inserts the read or write statements that are consistent with
c   the common blocks whose variables are intended to be read or written,
c   and thus create the 2 i/o subroutines.
c
c  This is intended to be a self-contained Fortran-77 program.
c
c  To run, simply compile with generic Fortran-77 compiler, and
c  execute.  No input or command lines are necessary.  E.g.
c  on unix:  f77 -o chekres chekres.f ; chekres 
c
c  Assumptions:
c    Common blocks are defined in the file <blokfil>.
c    Restart input read statements in the file <inresfil>.
c    Restart output write statements in the file <ouresfil>.
c    Common, read, write stmt keywords do not have imbedded blanks.
c    Variable names do not have imbedded blanks.
c    Variable names do not cross continuation line boundaries.
c    Variable & common block names are MAXNAME or less characters.
c    All characters will be mapped to lower case before analysis.
c    The character '!' signals the beginning of an in-line comment.
c    The common blocks named <noblok()> are not read/written in restart file.
c    "common /name/" is all on the same line.
c    "read(lun,...)" is all on the same line, & with no space before (.
c    "write(lun)" is all on the same line, & with no space before (.
c    The /aer1/ common block contains fundamental control parameters
c     such as itime, time, etc.
c
c  Expected reports from this program:
c    If no problems are found, the message: "All ok".
c    If problems, an attempt is made to describe the problem.
c
c  Things that control this program:
c    MAXBLOK, MAXVAR, NNOBLOK,
c    READSTR, WRITESTR
c    INSERT,
c    blokfil, intemplate, inresfil, outemplate, ouresfil,
c    noblok(),
c    
c
c
c  Define basic symbolic constants
c
      parameter(MAXBLOK=20)   ! max # common blocks
      parameter(MAXVAR=300)   ! max # vars in each common block
      parameter(MAXNAME=16)   ! max # chars in each common block var name
      parameter(LUNI=1)       ! lun for input files
      parameter(LUNO=2)       ! lun for output files
      parameter(NNOBLOK=2)    ! # common blocks that do not participate in restart
c
c
c  Declare & define restart read statement string (use all lower case)
c
      character*(*) READSTR
      parameter(READSTR='read(lunires')
c
c
c  Declare & define restart write statement string (use all lower case)
c
      character*(*) WRITESTR
      parameter(WRITESTR='write(lunores')
c
c
c  Declare & define template file string marking where to insert code
c
      character*(*) INSERT
      parameter(INSERT='<<< Insert restart i/o code here >>>')
c
c
c  Declare variables to hold common block names, variable names, etc
c
      integer nblok                                ! # common blocks
      integer nvar(MAXBLOK)                        ! # vars in each common block
      integer nvarres(MAXBLOK)                     ! # vars in each read/write block
      character*(MAXNAME) bloknam(MAXBLOK)         ! common block names
      character*(MAXNAME) varnam(MAXVAR,MAXBLOK)   ! var names in each block
      character*(MAXNAME) resnam(MAXVAR,MAXBLOK)   ! restart read/write names
      character*(40) blokfil                       ! common block include file name
      character*(40) inresfil                      ! file name of restart reader 
      character*(40) intemplate                    ! file name of restart reader template
      character*(40) ouresfil                      ! file name of restart writer
      character*(40) outemplate                    ! file name of restart writer template
      character*(MAXNAME) noblok(NNOBLOK)          ! common block not dumped to restart
c
c
c  Report startup of this program
c
      write(*,*) 'chekres now running ...'
c
c
c  Define names of input files
c
      blokfil = 'globaer.h' 
      intemplate = 'initres.template'
      inresfil = 'initres.f'
      outemplate = 'outres.template'
      ouresfil = 'outres.f'
c
c
c  Define names of non-restart common blocks
c
      noblok(1) = 'aer0'
      noblok(2) = 'aer0s'
c
c
c  Input & store the common block names & variable names in each common block
c
      call inblok
     $  (blokfil, MAXBLOK,MAXVAR,MAXNAME, LUNI,
     $   nblok, nvar, bloknam, varnam)
c
c
c  Create restart reader source code
c
      call genread
     $  (intemplate, inresfil, MAXBLOK,MAXVAR,MAXNAME, LUNI, LUNO,
     $   nblok, nvar, bloknam, varnam, NNOBLOK, noblok, INSERT)
c
c
c  Create restart writer source code
c
      call genwrite
     $  (outemplate, ouresfil, MAXBLOK,MAXVAR,MAXNAME, LUNI, LUNO,
     $   nblok, nvar, bloknam, varnam, NNOBLOK, noblok, INSERT)
c
c
c  Input & check restart reader source code for expected variable names
c
      call chkres
     $  (inresfil,READSTR, MAXBLOK,MAXVAR,MAXNAME, LUNI,
     $   nblok, nvar, bloknam, varnam, NNOBLOK, noblok, nvarres, resnam)
c
c
c  Input & check restart writer source code for variable names
c
      call chkres
     $  (ouresfil,WRITESTR, MAXBLOK,MAXVAR,MAXNAME, LUNI,
     $   nblok, nvar, bloknam, varnam, NNOBLOK, noblok, nvarres, resnam)
c
c
c  Report all ok if control reaches here
c
      write(*,*) 'chekres--All ok'
c
c
c  Terminate normally
c
      stop
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine inblok
     $  (filnam, maxblok, maxvar, maxname, lun,
     $   nblok, nvar, bloknam, varnam)
c
c
c  @(#) inblok.f  McKie  Oct-1995
c  This routine scans the source code of the file containing global
c  common blocks, looks for common block statements, and stores
c  the names of each common block, the number of variables in each
c  common block, and the names of the variables in each common block.
c
c  Input:
c            filnam = name of input global
c           maxblok = max # common blocks
c            maxvar = max # variables in each common block
c           maxname = max # characters in each common block or variable name
c               lun = i/o unit number to use
c
c  Output:
c             nblok = # common blocks
c           nvar(j) = # variables in j-th common block
c        bloknam(j) = name of j-th common block
c       varnam(i,j) = i-th variable name in j-th common block
c
c
c  Declare subprogram args
c
      character*(*) filnam
      integer maxblok
      integer maxvar
      integer maxname
      integer lun
      integer nblok
      integer nvar(maxblok)
      character*(*) bloknam(maxblok)
      character*(*) varnam(maxvar,maxblok)
c
c
c  Define local symbolic constants
c
      parameter(MAXLINE=72)
c
c
c  Declare local variables
c
      logical incommon
      logical ok
      character*(MAXLINE) line
c
c
c  Compute number of chars in file name
c
      call dblank(filnam, nfilnam)
c
c
c  Report beginning of common block input
c
      write(*,*) 'chekres--Examine common blocks from file: ',
     $  filnam(1:nfilnam)
c
c
c  Initialize number of common blocks
c
      nblok = 0
c
c
c  Initialize # variables in each block
c
      do j=1,maxblok
       nvar(j) = 0
      enddo
c
c
c  Initialize variable names in each block
c
      do j=1,maxblok
       do i=1,maxvar
        varnam(i,j) = ' '
       enddo
      enddo
c
c
c  Attempt to open the input file
c
      open(unit=lun,file=filnam,status='old',iostat=ios)
      if( ios .ne. 0 )then
       write(*,*) 'Error--Cant open input file ',filnam(1:nfilnam)
       stop
      endif
c
c
c  Initialize flag to indicate not within common block statement
c
      incommon = .false.
c
c
c  Initialize line counter
c
      kline = 0
c
c
c  Attempt to input next input file line
c
 2100 read(lun,'(a)',iostat=ios) line
c
c
c  Continue processing line unless eof was encountered
c
      if( ios .eq. 0 )then
c
c
c  Increment line counter
c
       kline = kline + 1
c
c
c  Find last non-blank char in line
c
       call dblank(line, nline)
c
c
c  If there is a Fortran internal comment, remove comment from line processing
c
       i = index(line,'!')
       if( i .gt. 0 ) nline = i - 1
c
c
c  Map chars in line to lower case
c
       do i=1,nline
        call tolower( line(i:i) )
       enddo
c
c
c  Ignore current line if it is an entire comment line
c
       if( line(1:1) .eq. 'c' ) goto 2100
c
c
c  Check to see if common stmt is beginning
c
       i0 = index( line(1:nline), 'common' )
       if( i0 .gt. 0 )then
        if( line(1:i0-1) .eq. ' ' )then
         incommon = .true.
c        write(*,*) 'Debug: ----------'
c        write(*,*) 'Debug: ' // line(1:nline)
         i1 = index( line(i0:nline), '/' )
         if( i1 .eq. 0 )then
          write(*,*) 'Error--At line ',kline,' of ',filnam(1:nfilnam)
          write(*,*) '       Expected common block /name/ missing'
          stop
         endif
         i1 = i0 + i1
         i2 = index( line(i1:nline), '/' )
         if( i2 .eq. 0 )then
          write(*,*) 'Error--At line ',kline,' of ',filnam(1:nfilnam)
          write(*,*) '       Expected common block /name/ missing'
          stop
         endif
         i2 = i1 + i2 - 2 
         if( ( i2 - i1 ) .lt. 0 )then
          write(*,*) 'Error--At line ',kline,' of ',filnam(1:nfilnam)
          write(*,*) '       Expected common block /name/ missing'
          stop
         endif
         nblok = nblok + 1
         if( nblok .gt. maxblok )then
          write(*,*) 'Error--Too many common blocks. Increase MAXBLOK'
          stop
         endif
         bloknam(nblok) = line(i1:i2)
         call rmblank( bloknam(nblok), ns)
         if( ns .lt. 1 )then
          write(*,*) 'Error--At line ',kline,' of ',filnam(1:nfilnam)
          write(*,*) '       Expected common block /name/ missing'
          stop
         endif
c        write(*,*) 'Debug: nblok=',nblok,'  bloknam=',bloknam(nblok)
         nvar(nblok) = 0
         m1 = i2 + 2
         m2 = nline
         ns = m2 - m1 + 1
         call getnam
     $     (maxvar,nvar(nblok),varnam(1,nblok),line(m1:m2),ns,ok)
         if( .not. ok )then
          write(*,*) 'Error--At line ',kline,' of ',filnam(1:nfilnam)
          write(*,*) '       Trouble recognizing variable name'
          stop
         endif
        endif
        goto 2100
       endif
c
c
c  If within common block stmt, check for continuation line
c
       if( incommon )then
        if( line(6:6) .eq. ' ' )then
         incommon = .false.
        else
c        write(*,*) 'Debug: ' // line(1:nline)
         m1 = 7
         m2 = nline
         ns = m2 - m1 + 1
         call getnam
     $     (maxvar,nvar(nblok),varnam(1,nblok),line(m1:m2),ns,ok)
         if( .not. ok )then
          write(*,*) 'Error--At line ',kline,' of ',filnam(1:nfilnam)
          write(*,*) '       Trouble recognizing variable name'
          stop
         endif
        endif
       endif
c
c
c  Go do next input line
c
       goto 2100
      endif
c
c
c  Close input file
c
      close(unit=lun)
c
c
c  Return to caller with list of common blocks & their variable names
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine genread
     $  (template, filnam, maxblok,maxvar,maxname, luni, luno,
     $   nblok, nvar, bloknam, varnam, nnoblok, noblok, insert )
c
c
c  @(#) genread.f  McKie  Jul-1996
c  This routine generates the restart file reader subroutine source code,
c  copying lines from a template file, adding code to do actual restart
c  file input, using one read statement for each common block that
c  is in the bloknam() list but not in the nnoblok() list.
c
c
c  Input:
c          template = name of input template file for restart input source code
c            filnam = name of output reader subroutine source file
c           maxblok = max # common blocks
c            maxvar = max # variables in each common block
c           maxname = max # characters in each common block or variable name
c              luni = input i/o unit number to use
c              luno = output i/o unit number to use
c             nblok = # common blocks
c           nvar(j) = # variables in j-th common block
c        bloknam(j) = name of j-th common block
c           nnoblok = # common blocks that are not read/written in restart
c         noblok(m) = name of m-th common block that is not read/written in restart
c       varnam(i,j) = i-th variable name in j-th common block
c            insert = inserted code goes after line with this string in template file
c
c
c  Declare subprogram args
c
      character*(*) template
      character*(*) filnam
      integer maxblok
      integer maxvar
      integer maxname
      integer luni
      integer luno
      integer nblok
      integer nvar(maxblok)
      character*(*) bloknam(maxblok)
      character*(*) varnam(maxvar,maxblok)
      integer nnoblok
      character*(*) noblok(nnoblok)
      character*(*) insert
c
c
c  Define local symbolic constants
c
      parameter(MAXLINE=90)
      parameter(MAXF77=72)
      parameter(NINDENT=7)
c
c
c  Declare local variables
c
      logical doneinsert
      logical useblok
      logical marker
      character*(MAXLINE) line
      character*(24) lastname
      character*(NINDENT) indent
      character*(MAXLINE) text
c
c
c  Define formats
c
    1 format(a)
    2 format('chekres--Error: Cant open template file ',a)
    3 format('chekres--Error: Cant find insert marker line in',
     $  ' template file',/,
     $  ' Template file: ',a,/,
     $  ' Marker string: ',a)
    4 format('chekres--Generating source file ',a,' from ',a)
c
c
c  Report source code file being generated
c
      call dblank(filnam, ns1)
      call dblank(template, ns2)
      write(*,4) filnam(1:ns1), template(1:ns2)
c
c
c  Define the <indent> string to be all blanks
c
      indent = ' '
c
c
c  Save copy of current reader subroutine source if it exists to backup file 
c
      open(unit=luni,file=filnam,status='old',err=2190)
       call dblank(filnam, ns)
       text = filnam(1:ns) // '.save'
       call dblank(text, ns)
       open(unit=luno,file=text(1:ns),status='unknown')
 2110   read(luni,1,iostat=ios) line
        if( ios .eq. 0 )then
         call dblank(line, nline)
         write(luno,1) line(1:nline)
         goto 2110
        endif
       close(unit=luni)
       close(unit=luno)
 2190 continue
c
c
c  Open input template file
c
      open(unit=luni,file=template,status='old',err=8110)
c
c
c  Open output reader subroutine file
c
      open(unit=luno,file=filnam,status='unknown')
c
c
c  Initialize flag to indicate that inserted code is not yet done
c
      doneinsert = .false.
c
c
c  Input template lines until eof
c
 4100 read(luni,1,end=7100) line
c
c
c  Compute number of chars to last non-blank char in current template line
c
       call dblank(line, nline)
c
c
c  Define flag for if current template line is marker line to insert i/o stmts
c
       marker = index(line,insert) .gt. 0
c
c
c  Copy current template line to output subroutine source file unless it is marker
c
       if( .not. marker ) write(luno,1) line(1:nline)
c
c
c  Insert read stmts into source if current template line is place to do this
c
       if( marker )then
c
c
c  Indicate that code was inserted into template
c
        doneinsert = .true.
c
c
c  Visit each of the known common blocks
c
        do iblok=1,nblok
c
c
c  Check if current common block name is in list of those not in restart file
c
         useblok = .true.
         do k=1,nnoblok
          if( bloknam(iblok) .eq. noblok(k) ) useblok = .false.
         enddo
c
c
c  Add i/o stmts for current common block if it is a restart common block
c
         if( useblok )then
c
c
c  Add the read stmt containing all variables in current common block
c
          call dblank(bloknam(iblok), ns)
c         write(*,*) 'Debug. Add read code for ', bloknam(iblok)(1:ns)
          write(LUNO,1) 'c'
          write(LUNO,1) 'c'
          text = 'c  Input contents of common /' //
     $     bloknam(iblok)(1:ns) // '/ & check for safety marker value'
          call dblank(text, ns)
          write(LUNO,1) text(1:ns)
          write(LUNO,1) 'c'
          write(LUNO,1) indent // 'read(LUNIRES,end=8100,err=8100)'
          m = MAXF77 + 1
          nvarlin = 0
          do ivar=1,nvar(iblok)
           call dblank(varnam(ivar,iblok), ns)
           mnew = m + 1 + ns
           if( ivar .lt. nvar(iblok) ) mnew = mnew + 1
           if( mnew .gt. MAXF77 )then
            if( nvarlin .gt. 0 )then
             call dblank(line, nline)
             write(LUNO,1) line(1:nline)
            endif
            line = '     $'
            m = NINDENT + 2 
           endif
           line(m+1:) = ' ' // varnam(ivar,iblok)(1:ns)
           nvarlin = nvarlin + 1
           m = m + 1 + ns
           if( ivar .lt. nvar(iblok) )then
            m = m + 1
            line(m:m) = ','
           endif
c          write(*,*) 'Debug. ivar, m, line=',ivar,m,line(1:m)
          enddo
          if( nvarlin .gt. 0 )then
           call dblank(line, nline)
           write(LUNO,1) line(1:nline)
          endif
c
c
c  Add code to check value of last special variable in common block
c   (Use different code for numeric or character common blocks)
c
          lastname = varnam(nvar(iblok),iblok)
          call dblank(lastname, nl)
          call dblank(bloknam(iblok),ns)
          if( bloknam(iblok)(ns:ns) .eq. 's' )then
           write(LUNO,1) indent //
     $       'if( ' // lastname(1:nl) // ' .ne. CSAFETY )then'
           write(LUNO,1) indent //
     $       ' write(LUNOPRT,4) ''' // lastname(1:nl) // ''', ' //
     $       lastname(1:nl)
           write(LUNO,1) indent //
     $       ' write(LUNOPRT,10)'
           write(LUNO,1) indent //
     $       ' stop'
           write(LUNO,1) indent //
     $       'endif'
          else
           write(LUNO,1) indent //
     $       'if( ' // lastname(1:nl) // ' .ne. ISAFETY )then'
           write(LUNO,1) indent //
     $       ' write(LUNOPRT,3) ''' // lastname(1:nl) // ''', ' //
     $       lastname(1:nl)
           write(LUNO,1) indent //
     $       ' write(LUNOPRT,10)'
           write(LUNO,1) indent //
     $       ' stop'
           write(LUNO,1) indent //
     $       'endif'
          endif
c
c
c  End of <useblok> conditional
c
         endif
c
c
c  Go do next known common block
c
        enddo
c
c
c  End of inserted code
c
       endif
c
c
c  Go input next template line
c
       goto 4100
 7100 continue
c
c
c  Close input template & output reader subroutine files
c
      close(unit=luni)
      close(unit=luno)
c
c
c  Ensure that code was inserted into template file
c
      if( .not. doneinsert )then 
       call dblank(template, ns)
       write(*,3) template(1:ns), insert
       stop
      endif
c
c
c  Return to caller with restart reader subroutine source generated
c
      return
c
c
c  Error processing
c
 8110 continue
      call dblank(template, ns)
      write(*,2) template(1:ns)
      stop
c
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine genwrite
     $  (template, filnam, maxblok,maxvar,maxname, luni, luno,
     $   nblok, nvar, bloknam, varnam, nnoblok, noblok, insert )
c
c
c  @(#) genwrite.f  McKie  Jul-1996
c  This routine generates the restart file writer subroutine source code,
c  copying lines from a template file, adding code to do actual restart
c  file output, using one write statement for each common block that
c  is in the bloknam() list but not in the nnoblok() list.
c
c
c  Input:
c          template = name of input template file for restart output source code
c            filnam = name of output writer subroutine source file
c           maxblok = max # common blocks
c            maxvar = max # variables in each common block
c           maxname = max # characters in each common block or variable name
c              luni = input i/o unit number to use
c              luno = output i/o unit number to use
c             nblok = # common blocks
c           nvar(j) = # variables in j-th common block
c        bloknam(j) = name of j-th common block
c           nnoblok = # common blocks that are not read/written in restart
c         noblok(m) = name of m-th common block that is not read/written in restart
c       varnam(i,j) = i-th variable name in j-th common block
c            insert = inserted code goes after line with this string in template file
c
c
c  Declare subprogram args
c
      character*(*) template
      character*(*) filnam
      integer maxblok
      integer maxvar
      integer maxname
      integer luni
      integer luno
      integer nblok
      integer nvar(maxblok)
      character*(*) bloknam(maxblok)
      character*(*) varnam(maxvar,maxblok)
      integer nnoblok
      character*(*) noblok(nnoblok)
      character*(*) insert
c
c
c  Define local symbolic constants
c
      parameter(MAXLINE=90)
      parameter(MAXF77=72)
      parameter(NINDENT=7)
c
c
c  Declare local variables
c
      logical doneinsert
      logical useblok
      logical marker
      character*(MAXLINE) line
      character*(24) lastname
      character*(NINDENT) indent
      character*(MAXLINE) text
c
c
c  Define formats
c
    1 format(a)
    2 format('chekres--Error: Cant open template file ',a)
    3 format('chekres--Error: Cant find insert marker line in',
     $  ' template file',/,
     $  ' Template file: ',a,/,
     $  ' Marker string: ',a)
    4 format('chekres--Generating source file ',a,' from ',a)
c
c
c  Report source code file being generated
c
      call dblank(filnam, ns1)
      call dblank(template, ns2)
      write(*,4) filnam(1:ns1), template(1:ns2)
c
c
c  Define the <indent> string to be all blanks
c
      indent = ' '
c
c
c  Save copy of current writer subroutine source if it exists to backup file 
c
      open(unit=luni,file=filnam,status='old',err=2190)
       call dblank(filnam, ns)
       text = filnam(1:ns) // '.save'
       call dblank(text, ns)
       open(unit=luno,file=text(1:ns),status='unknown')
 2110   read(luni,1,iostat=ios) line
        if( ios .eq. 0 )then
         call dblank(line, nline)
         write(luno,1) line(1:nline)
         goto 2110
        endif
       close(unit=luni)
       close(unit=luno)
 2190 continue
c
c
c  Open input template file
c
      open(unit=luni,file=template,status='old',err=8110)
c
c
c  Open output writer subroutine file
c
      open(unit=luno,file=filnam,status='unknown')
c
c
c  Initialize flag to indicate that inserted code is not yet done
c
      doneinsert = .false.
c
c
c  Input template lines until eof
c
 4100 read(luni,1,end=7100) line
c
c
c  Compute number of chars to last non-blank char in current template line
c
       call dblank(line, nline)
c
c
c  Define flag for if current template line is marker line to insert i/o stmts
c
       marker = index(line,insert) .gt. 0
c
c
c  Copy current template line to output subroutine source file unless it is marker
c
       if( .not. marker ) write(luno,1) line(1:nline)
c
c
c  Insert write stmts into source if current template line is place to do this
c
       if( marker )then
c
c
c  Indicate that code was inserted into template
c
        doneinsert = .true.
c
c
c  Visit each of the known common blocks
c
        do iblok=1,nblok
c
c
c  Check if current common block name is in list of those not in restart file
c
         useblok = .true.
         do k=1,nnoblok
          if( bloknam(iblok) .eq. noblok(k) ) useblok = .false.
         enddo
c
c
c  Add i/o stmts for current common block if it is a restart common block
c
         if( useblok )then
c
c
c  Add the comment stmts to describe output of current common block
c
          call dblank(bloknam(iblok), ns)
c         write(*,*) 'Debug. Add write code for ', bloknam(iblok)(1:ns)
          write(LUNO,1) 'c'
          write(LUNO,1) 'c'
          text = 'c  Load safety value' //
     $     ' & output contents of common /' //
     $     bloknam(iblok)(1:ns) // '/'
          call dblank(text, ns)
          write(LUNO,1) text(1:ns)
          write(LUNO,1) 'c'
c
c
c  Add code to load value of last special variable in common block
c   (Use different code for numeric or character common blocks)
c
          lastname = varnam(nvar(iblok),iblok)
          call dblank(lastname, nl)
          call dblank(bloknam(iblok),ns)
          if( bloknam(iblok)(ns:ns) .eq. 's' )then
           write(LUNO,1) indent //
     $       lastname(1:nl) // ' = CSAFETY'
          else
           write(LUNO,1) indent //
     $       lastname(1:nl) // ' = ISAFETY'
          endif
c
c
c  Add the write stmt containing all variables in current common block
c
          call dblank(bloknam(iblok), ns)
          write(LUNO,1) indent // 'write(LUNORES)'
          m = MAXF77 + 1
          nvarlin = 0
          do ivar=1,nvar(iblok)
           call dblank(varnam(ivar,iblok), ns)
           mnew = m + 1 + ns
           if( ivar .lt. nvar(iblok) ) mnew = mnew + 1
           if( mnew .gt. MAXF77 )then
            if( nvarlin .gt. 0 )then
             call dblank(line, nline)
             write(LUNO,1) line(1:nline)
            endif
            line = '     $'
            m = NINDENT + 2 
           endif
           line(m+1:) = ' ' // varnam(ivar,iblok)(1:ns)
           nvarlin = nvarlin + 1
           m = m + 1 + ns
           if( ivar .lt. nvar(iblok) )then
            m = m + 1
            line(m:m) = ','
           endif
c          write(*,*) 'Debug. ivar, m, line=',ivar,m,line(1:m)
          enddo
          if( nvarlin .gt. 0 )then
           call dblank(line, nline)
           write(LUNO,1) line(1:nline)
          endif
c
c
c  End of <useblok> conditional
c
         endif
c
c
c  Go do next known common block
c
        enddo
c
c
c  End of inserted code
c
       endif
c
c
c  Go input next template line
c
       goto 4100
 7100 continue
c
c
c  Close input template & output writer subroutine files
c
      close(unit=luni)
      close(unit=luno)
c
c
c  Ensure that code was inserted into template file
c
      if( .not. doneinsert )then 
       call dblank(template, ns)
       write(*,3) template(1:ns), insert
       stop
      endif
c
c
c  Return to caller with restart writer subroutine source generated
c
      return
c
c
c  Error processing
c
 8110 continue
      call dblank(template, ns)
      write(*,2) template(1:ns)
      stop
c
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine chkres
     $  (filnam, stmt, maxblok,maxvar,maxname, lun,
     $   nblok, nvar, bloknam, varnam, nnoblok, noblok, nvarres, resnam)
c
c
c  @(#) chkres.f  McKie  Oct-1995
c  This routine scans the source code of the file containing 
c  restart file input code or restart file output code,
c  and builds a list of input or output blocks of variables
c  that are assumed to be the equivalent of global common blocks.
c  It then compares this list with the list of known
c  common blocks and the variables in each block, and ensures that
c  the input or output blocks are identical to the known common blocks,
c  and that all appropriate common blocks are represented in the input
c  or output blocks of variables.
c
c
c  Input:
c            filnam = name of input global
c              stmt = 'read' for restart input, 'write' for restart output
c           maxblok = max # common blocks
c            maxvar = max # variables in each common block
c           maxname = max # characters in each common block or variable name
c               lun = i/o unit number to use
c             nblok = # common blocks
c           nvar(j) = # variables in j-th common block
c        bloknam(j) = name of j-th common block
c           nnoblok = # common blocks that are not read/written in restart
c         noblok(m) = name of m-th common block that is not read/written in restart
c       varnam(i,j) = i-th variable name in j-th common block
c        nvarres(j) = # variables in j-th common block
c       resnam(i,j) = array to hold i-th variable name in j-th read or write stmt
c
c
c  Declare subprogram args
c
      character*(*) filnam
      character*(*) stmt
      integer maxblok
      integer maxvar
      integer maxname
      integer lun
      integer nblok
      integer nvar(maxblok)
      character*(*) bloknam(maxblok)
      character*(*) varnam(maxvar,maxblok)
      integer nnoblok
      character*(*) noblok(nnoblok)
      integer nvarres(maxblok)
      character*(*) resnam(maxvar,maxblok)
c
c
c  Define local symbolic constants
c
      parameter(MAXLINE=72)
c
c
c  Declare local variables
c
      logical inrw
      logical ok
      character*(MAXLINE) line
c
c
c  Compute number of chars in file name
c
      call dblank(filnam, nfilnam)
c
c
c  Report beginning of common block input
c
      write(*,*) 'chekres--Checking source file: ', filnam(1:nfilnam)
c
c
c  Initialize number of read/write blocks
c
      nrwb = 0
c
c
c  Initialize # variables in each read/write block
c
      do j=1,maxblok
       nvarres(j) = 0
      enddo
c
c
c  Initialize variable names in each read/write block
c
      do j=1,maxblok
       do i=1,maxvar
        resnam(i,j) = ' '
       enddo
      enddo
c
c
c  Attempt to open the input file
c
      open(unit=lun,file=filnam,status='old',iostat=ios)
      if( ios .ne. 0 )then
       write(*,*) 'Error--Cant open input file ',filnam(1:nfilnam)
       stop
      endif
c
c
c  Initialize flag to indicate not within read/write statement
c
      inrw = .false.
c
c
c  Initialize line counter
c
      kline = 0
c
c
c  Attempt to input next input file line
c
 2100 read(lun,'(a)',iostat=ios) line
c
c
c  Continue processing line unless eof was encountered
c
      if( ios .eq. 0 )then
c
c
c  Increment line counter
c
       kline = kline + 1
c
c
c  Find last non-blank char in line
c
       call dblank(line, nline)
c
c
c  If there is a Fortran internal comment, remove comment from line processing
c
       i = index(line,'!')
       if( i .gt. 0 ) nline = i - 1
c
c
c  Map chars in line to lower case
c
       do i=1,nline
        call tolower( line(i:i) )
       enddo
c
c
c  Ignore current line if it is an entire comment line
c
       if( line(1:1) .eq. 'c' ) goto 2100
c
c
c  Check to see if read/write stmt is beginning
c
       i0 = index( line(1:nline), stmt )
       if( i0 .gt. 0 )then
        if( ( i0 .eq. 7 ) .or. ( line(7:i0-1) .eq. ' ' ) )then
         inrw = .true.
c        write(*,*) 'Debug: ----------'
c        write(*,*) 'Debug: ' // line(1:nline)
         i2 = index( line(i0:nline), ')' )
         if( i2 .eq. 0 )then
          write(*,*) 'Error--At line ',kline,' of ',filnam(1:nfilnam)
          write(*,*) '       Expected ', stmt, ' (lun) missing'
          stop
         endif
         i2 = i0 + i2 - 1 
         nrwb = nrwb + 1
         if( nrwb .gt. maxblok )then
          write(*,*) 'Error--Too many ', stmt,
     $               ' blocks. Increase MAXBLOK'
          stop
         endif
         nvarres(nrwb) = 0
         m1 = i2 + 1
         m2 = nline
         ns = m2 - m1 + 1
         call getnam
     $     (maxvar,nvarres(nrwb),resnam(1,nrwb),line(m1:m2),ns,ok)
         if( .not. ok )then
          write(*,*) 'Error--At line ',kline,' of ',filnam(1:nfilnam)
          write(*,*) '       Trouble recognizing variable name'
          stop
         endif
        endif
        goto 2100
       endif
c
c
c  If within read/write stmt, check for continuation line
c
       if( inrw )then
        if( line(6:6) .eq. ' ' )then
         inrw = .false.
        else
c        write(*,*) 'Debug: ' // line(1:nline)
         m1 = 7
         m2 = nline
         ns = m2 - m1 + 1
         call getnam
     $     (maxvar,nvarres(nrwb),resnam(1,nrwb),line(m1:m2),ns,ok)
         if( .not. ok )then
          write(*,*) 'Error--At line ',kline,' of ',filnam(1:nfilnam)
          write(*,*) '       Trouble recognizing variable name'
          stop
         endif
        endif
       endif
c
c
c  Go do next input line
c
       goto 2100
      endif
c
c
c  Close input file
c
      close(unit=lun)
c
c
c  Report problem is not same number of restart common blocks as
c  there are restart read/write blocks
c
      if( ( nblok - nnoblok ) .ne. nrwb )then
       write(*,*) 'chekres--Problem in restart ', stmt,
     $            ' source file ', filnam(1:nfilnam) 
       write(*,*) '  ', nblok - nnoblok, ' restart common blocks'
       write(*,*) '  ', nrwb, ' restart ', stmt, ' blocks'
       stop
      endif
c
c
c  Initialize counter of read/write blocks
c
      n = 0
c
c
c  Scan list of expected common blocks
c
      do 7100 j=1,nblok
c
c
c  Ignore current common block if it does not participate in restart
c
       do m=1,nnoblok
        if( bloknam(j) .eq. noblok(m) ) goto 7100
       enddo
c
c
c  Point to next expected read/write block index
c
       n = n + 1
c
c
c  Report problem if # variables in common block not same as read/write block
c
        if( nvar(j) .ne. nvarres(n) )then
         write(*,*) 'chekres--Problem in restart ', stmt,
     $              ' source file ', filnam(1:nfilnam) 
         write(*,*) '  ', nvar(j), ' variables in common block ',
     $              bloknam(j)
         write(*,*) '  ', nvarres(n), ' variables in ', stmt,
     $              ' stmt #', n
c        write(*,*) '  ', stmt, ' stmt should have variables:'
c        do i=1,nvar(j)
c         write(*,*) '    ', varnam(i,j)
c        enddo
         stop
c
c
c  Report if variable name list is not identical in common & read/write blocks
c
        else
         do i=1,nvar(j)
          if( varnam(i,j) .ne. resnam(i,n) )then
           write(*,*) 'Problem in restart ', stmt,
     $                ' source file ', filnam(1:nfilnam) 
           write(*,*) '  Variable names do not match ',
     $                ' for common block ', bloknam(j)
           write(*,*) '  and ', stmt, ' stmt # ',n
           do m=1,nvar(j)
            call dblank( varnam(i,j), ns1)
            call dblank( resnam(i,n), ns2)
            write(*,*) 'Variable # ',m,
     $                 '  common=',varnam(i,j)(1:ns1),
     $                 '  ', stmt, '=',resnam(i,n)(1:ns2)
           enddo
           stop
          endif
         enddo
        endif
c
c
c  Go examine next common block
c
 7100 continue
c
c
c  Return to caller with restart input or output checked
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine getnam(maxvar,n,vars,s,ns,ok)
c
c
c  @(#) getnam.f  McKie  Oct-1995
c  This routine appends the variable names from a string to a
c  list of variable names.
c
c  Input:
c            maxvar = max # names that can be in vars()
c                 n = current # names in vars()
c            vars() = current list of names
c                 s = string containing new names
c                ns = # chars in s
c
c  Output:
c                 n = new # names in vars()
c            vars() = new list of names
c                ok = .true. if no problems in appending more names
c
c
c  Declare subprogram args
c
      integer maxvar
      integer n
      character*(*) vars(*)
      character*(*) s
      integer ns
      logical ok
c
c
c  Declare local variables
c
      logical inname
      logical inparen
      logical isnamec
      logical isalpha
      character*(1) c
c
c
c  Initialize flag to indicate that we're not currently within a name
c
      inname = .false.
c
c
c  Initialize flag to indicate that we're not currently within paretheses
c
      inparen = .false.
c
c
c  Scan the chars in the <s> string
c
      do i=1,ns
c
c
c  Define current character within the <s> string
c
       c = s(i:i)
c
c
c  Check if current character is a variable name character
c
       call namchar(c, isnamec)
c
c
c  Check if current character is an alphabetic character
c
       call alpha(c, isalpha)
c
c
c  Start new variable name if not in parentheses, not in name, &
c  current char is an alphabetic
c
       if( (.not.inparen) .and. (.not.inname) .and. (isalpha) )then
        i1 = i
        inname = .true.
       endif
c
c
c  Add current char to name if not in parentheses, in name, & char is name char
c
       if( (.not.inparen) .and. (inname) .and. (isnamec) )then
        i2 = i
       endif
c
c
c  End variable if not in parentheses, in name, and either
c  char is not a name char, or char is a name char & last char on line
c
       if( (.not.inparen) .and. (inname) )then
        if( (.not.isnamec) .or. ( (isnamec) .and. (i.eq.ns) ) )then
         inname = .false.
         n = n + 1
         if( n .gt. maxvar )then
          write(*,*) 'Error--Too many names. Increase MAXVAR'
          ok = .false.
          return
         else
          vars(n) = s(i1:i2)
          call rmblank( vars(n), nvars )
c         write(*,*) 'Debug: n=', n, ' New name=', vars(n)(1:nvars)
         endif
        endif
       endif
c
c
c  Begin parentheses if char is '('
c
       if( c .eq. '(' )then
        inparen = .true.
       endif
c
c
c  End parentheses if char is ')'
c
       if( c .eq. ')' )then
        inparen = .false.
       endif
c
c
c  Go process next character in s string
c
      enddo
c
c
c  Set result flag to indicate no problems in appending variable names
c
      ok = .true.
c
c
c  Return to caller with variable names from <s> appended to list in <vars()>
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine namchar(c, isnamec)
c
c
c  @(#) namchar.f  McKie  Oct-1995
c  This routine checks a character to see if it is a character that
c  could be within a Fortran name or not.
c
c
c  Input:
c                 c = character to be examined
c
c  Output:
c           isnamec = .true. if <c> could be a char in a f77 name, else .false.
c
c
c  Declare subprogram args
c
      character*(1) c
      logical isnamec
c
c
c  Compute whether char could be a char in a variable name or not
c
      isnamec =
     $    ( ( c .ge. 'a' ) .and. ( c .le. 'z' ) ) .or.
     $    ( ( c .ge. '0' ) .and. ( c .le. '9' ) ) .or.
     $    ( ( c .ge. 'A' ) .and. ( c .le. 'Z' ) ) .or.
     $    ( c .eq. '_' ) 
c
c
c  Return to caller with character examined & <isnamec> set
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine alpha(c, isalpha)
c
c
c  @(#) alpha.f  McKie  Oct-1995
c  This routine checks a character to see if it is an alphabetic char.
c
c
c  Input:
c                 c = character to be examined
c
c  Output:
c           isalpha = .true. if <c> is an alphabetic char, else .false.
c
c
c  Declare subprogram args
c
      character*(1) c
      logical isalpha
c
c
c  Compute if char is an alphabetic or not
c
      isalpha =
     $    ( ( c .ge. 'a' ) .and. ( c .le. 'z' ) ) .or.
     $    ( ( c .ge. 'A' ) .and. ( c .le. 'Z' ) ) 
c
c
c  Return to caller with character examined & <isalpha> set
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine tolower(c)
c
c
c  @(#) tolower.f  McKie  Oct-1995
c  This routine examines an ascii character, and if it is an upper
c  case alpha character, it maps it to the lower case version of that
c  character.
c
c  Input:
c           c = character to be examined
c
c  Output:
c           c = original character if not upper case, else lower case equiv
c
c
c  Declare subprogram args
c
      character*(1) c
c
c
      if( ( c .ge. 'A' ) .and. ( c .le. 'Z' ) )then
       ic = ichar(c)
       ic6 = ic / 32
       if( ( ( ic6 / 2 ) * 2 ) .eq. ic6 )then
        c = char( ic + 32 )
       endif
      endif
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine dblank(s, ns)
c
c
c  @(#) dblank.f  McKie  Jun-1988
c  This routine finds index of last nonblank char in a string.
c
c  Argument list input:
c    s = string to be examined
c
c  Argument list output:
c    ns = # chars in s, up to & including last non-blank
c
c
c  Declare subprogram arg(s)
c
      character*(*) s
      integer ns
c
c
c  Find last non-blank char in string (beyond 1st char), return to caller
c
      ns = len(s)
 2100 if( ( s(ns:ns) .eq. ' ' ) .and. ( ns .gt. 1 ) )then
       ns = ns - 1
       goto 2100
      endif
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine rmblank(s, ns)
c
c
c  @(#) rmblank.f  McKie  Jul-1988
c  This routine removes all blanks from a string & left justifies it.
c
c  Input:    s = string to be examined
c
c  Output:
c           ns = # non-blank chars in s
c            s = string with all blanks removed, left justified
c
c
c  Declare subprogram arg(s)
c
      character*(*) s
c
c
c  Declare local variables
c
      character*1 c
c
c
c  Move non-blank chars to left of word, squeezing out blanks
c
      nsold = len(s)
      ns = 0
      do 2100 i=1,nsold
       c = s(i:i)
       if( c .ne. ' ' )then
        ns = ns + 1
        s(ns:ns) = c
       endif
 2100 continue
c
c
c  Fill out remaining rightmost chars in string with blanks
c
      do 2300 i=ns+1,nsold
       s(i:i) = ' '
 2300 continue
c
c
c  Return to caller with blanks removed, & non-blank chars left justified
c
      return
      end

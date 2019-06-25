c  Test program
      open(unit=10,file='test.out',status='unknown',
     $     form='unformatted')
      write(10) NZ
      close(10)
      stop
      end

      function cvmgt(x1,x2,x3)
c
c   Include implicit declarations
c
      include 'precision.h'
 
      logical x3

      if (x3) then
        cvmgt = x1
      else
        cvmgt = x2
      endif

      return
      end

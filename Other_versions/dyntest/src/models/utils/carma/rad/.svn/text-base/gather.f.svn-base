      subroutine gather(n,a,b,index)
c
c   Include implicit declarations
c
      include 'precision.h'
 
      integer n,index(n)
      dimension a(n),b(n)

      do 10 i = 1, n
        a(i) = b(index(i))
   10 continue

      return
      end

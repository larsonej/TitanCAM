	  �6  �   k820309    �          12.0        &mBO                                                                                                           
       /home/larsonej/titancam/src/models/atm/cam/src/physics/cam1/volcanicmass.F90 VOLCANICMASS              READ_VOLCANIC_MASS GET_VOLCANIC_MASS VOLCANIC_INITIALIZE                  � @                              
       R8 SHR_KIND_R8                  � @                               
       PLAT PLON PLOND MASTERPROC PLEV PLEVP                                                     
       LATDEG          @       � @                               
       BEGCHUNK ENDCHUNK PCOLS PVER PVERP                      @                              
       SCATTER_FIELD_TO_CHUNK GET_NCOLS_P                      @                              
       BNDTVVOLC                  � @ @                             
       ENDRUN                      @                              
       VALIDFACTORS                  � @                         	     
                     � @                          
     '                     #DAY    #TOD                � D                                                           �                                                         #ESMF_TOD                   � @                              '                    #TYPE    #SEC    #MSEC                 � D                                                            � D                                                           � D                                                            � @                               '                    #CALENDAR    #YEAR    #MONTH    #DAY    #TOD    #JULIANDAY    #DAYOFYEAR                �                                    �                      #ESMF_CALENDAR                   � @                              '�                    #TYPE    #DIM    #DIMRUNNINGSUM    #DIY                 � D                                                            � D                                                           p          p            p                                       � D                                        p                   p          p            p                                       � D                                 �                         � D                                 �                         � D                                 �                         � D                                 �                         �                                           �              #ESMF_TOD                � D                                                         � D                                                           � @                               '�                   #NSTEP    #STEPSIZE     #STARTDATE !   #STOPDATE "   #BASEDATE #   #CURRDATE $   #PREVDATE %               � D                                                           �                                                           #ESMF_TIME 
               �                               !            (              #ESMF_DATE                �                               "            H             #ESMF_DATE                �                               #            h             #ESMF_DATE                �                               $            �             #ESMF_DATE                �                               %            �             #ESMF_DATE               @   !                           &            #         @                                 '                  #ENDRUN%PRESENT (   #MSG )                 @                            (     PRESENT           
  @                             )                    1 %         @                               *                          #VALIDFACTORS%ABS +   #FACT1 ,   #FACT2 -                 @                            +     ABS           
   @                             ,     
                
   @                             -     
      #         @                                  .                   #READ_VOLCANIC_MASS%ESMF_CALENDAR /   #READ_VOLCANIC_MASS%ESMF_DATE 4   #READ_VOLCANIC_MASS%ESMF_TOD :   #READ_VOLCANIC_MASS%ESMF_TIME @   #READ_VOLCANIC_MASS%ESMF_TIMEMGR C                  � @                         /     '�                    #TYPE 0   #DIM 1   #DIMRUNNINGSUM 2   #DIY 3                � D                            0                                � D                            1                               p          p            p                                       � D                            2            p                   p          p            p                                       � D                            3     �                           � @                         4     '                    #CALENDAR 5   #YEAR 6   #MONTH 7   #DAY 8   #TOD 9   #JULIANDAY >   #DAYOFYEAR ?               �                               5     �                      #READ_VOLCANIC_MASS%ESMF_CALENDAR /               � D                            6     �                         � D                            7     �                         � D                            8     �                         �                               9            �              #READ_VOLCANIC_MASS%ESMF_TOD :                  � @                         :     '                    #TYPE ;   #SEC <   #MSEC =                � D                            ;                                � D                            <                               � D                            =                              � D                            >                             � D                            ?                               � @                         @     '                     #DAY A   #TOD B               � D                            A                               �                               B                          #READ_VOLCANIC_MASS%ESMF_TOD :                 � @                          C     '�                   #NSTEP D   #STEPSIZE E   #STARTDATE F   #STOPDATE G   #BASEDATE H   #CURRDATE I   #PREVDATE J               � D                            D                               �                               E                           #READ_VOLCANIC_MASS%ESMF_TIME @               �                               F            (              #READ_VOLCANIC_MASS%ESMF_DATE 4               �                               G            H             #READ_VOLCANIC_MASS%ESMF_DATE 4               �                               H            h             #READ_VOLCANIC_MASS%ESMF_DATE 4               �                               I            �             #READ_VOLCANIC_MASS%ESMF_DATE 4               �                               J            �             #READ_VOLCANIC_MASS%ESMF_DATE 4   #         @                                  K                  #GET_VOLCANIC_MASS%ESMF_CALENDAR L   #GET_VOLCANIC_MASS%ESMF_DATE Q   #GET_VOLCANIC_MASS%ESMF_TOD W   #GET_VOLCANIC_MASS%ESMF_TIME ]   #GET_VOLCANIC_MASS%ESMF_TIMEMGR `   #LCHNK h   #AEROSOL i                  � @                         L     '�                    #TYPE M   #DIM N   #DIMRUNNINGSUM O   #DIY P                � D                            M                                � D                            N                               p          p            p                                       � D                            O            p                   p          p            p                                       � D                            P     �                           � @                         Q     '                    #CALENDAR R   #YEAR S   #MONTH T   #DAY U   #TOD V   #JULIANDAY [   #DAYOFYEAR \               �                               R     �                      #GET_VOLCANIC_MASS%ESMF_CALENDAR L               � D                            S     �                         � D                            T     �                         � D                            U     �                         �                               V            �              #GET_VOLCANIC_MASS%ESMF_TOD W                  � @                         W     '                    #TYPE X   #SEC Y   #MSEC Z                � D                            X                                � D                            Y                               � D                            Z                              � D                            [                             � D                            \                               � @                         ]     '                     #DAY ^   #TOD _               � D                            ^                               �                               _                          #GET_VOLCANIC_MASS%ESMF_TOD W                 � @                          `     '�                   #NSTEP a   #STEPSIZE b   #STARTDATE c   #STOPDATE d   #BASEDATE e   #CURRDATE f   #PREVDATE g               � D                            a                               �                               b                           #GET_VOLCANIC_MASS%ESMF_TIME ]               �                               c            (              #GET_VOLCANIC_MASS%ESMF_DATE Q               �                               d            H             #GET_VOLCANIC_MASS%ESMF_DATE Q               �                               e            h             #GET_VOLCANIC_MASS%ESMF_DATE Q               �                               f            �             #GET_VOLCANIC_MASS%ESMF_DATE Q               �                               g            �             #GET_VOLCANIC_MASS%ESMF_DATE Q             
  @@                              h                     
D  @                             i                   
               &                   &                                           #         @                                  j                   #VOLCANIC_INITIALIZE%ESMF_CALENDAR k   #VOLCANIC_INITIALIZE%ESMF_DATE p   #VOLCANIC_INITIALIZE%ESMF_TOD v   #VOLCANIC_INITIALIZE%ESMF_TIME |   #VOLCANIC_INITIALIZE%ESMF_TIMEMGR                   � @                         k     '�                    #TYPE l   #DIM m   #DIMRUNNINGSUM n   #DIY o                � D                            l                                � D                            m                               p          p            p                                       � D                            n            p                   p          p            p                                       � D                            o     �                           � @                         p     '                    #CALENDAR q   #YEAR r   #MONTH s   #DAY t   #TOD u   #JULIANDAY z   #DAYOFYEAR {               �                               q     �                      #VOLCANIC_INITIALIZE%ESMF_CALENDAR k               � D                            r     �                         � D                            s     �                         � D                            t     �                         �                               u            �              #VOLCANIC_INITIALIZE%ESMF_TOD v                  � @                         v     '                    #TYPE w   #SEC x   #MSEC y                � D                            w                                � D                            x                               � D                            y                              � D                            z                             � D                            {                               � @                         |     '                     #DAY }   #TOD ~               � D                            }                               �                               ~                          #VOLCANIC_INITIALIZE%ESMF_TOD v                 � @                               '�                   #NSTEP �   #STEPSIZE �   #STARTDATE �   #STOPDATE �   #BASEDATE �   #CURRDATE �   #PREVDATE �               � D                            �                               �                               �                           #VOLCANIC_INITIALIZE%ESMF_TIME |               �                               �            (              #VOLCANIC_INITIALIZE%ESMF_DATE p               �                               �            H             #VOLCANIC_INITIALIZE%ESMF_DATE p               �                               �            h             #VOLCANIC_INITIALIZE%ESMF_DATE p               �                               �            �             #VOLCANIC_INITIALIZE%ESMF_DATE p               �                               �            �             #VOLCANIC_INITIALIZE%ESMF_DATE p      �   b      fn#fn "     I   b   uapp(VOLCANICMASS    K  O   J  SHR_KIND_MOD    �  f   J  PMGRID       G   J  COMMAP    G  c   J  PPGRID    �  c   J  PHYS_GRID      J   J  FILENAMES    W  G   J  ABORTUTILS    �  M   J  TIMEINTERP    �  @   J  MPISHORTHAND '   +  b       ESMF_TIME+ESMF_TIMEMOD /   �  H   %   ESMF_TIME%DAY+ESMF_TIMEMOD=DAY +   �  ^   a   ESMF_TIME%TOD+ESMF_TIMEMOD %   3  m      ESMF_TOD+ESMF_TODMOD /   �  H   %   ESMF_TOD%TYPE+ESMF_TODMOD=TYPE -   �  H   %   ESMF_TOD%SEC+ESMF_TODMOD=SEC /   0  H   %   ESMF_TOD%MSEC+ESMF_TODMOD=MSEC '   x  �       ESMF_DATE+ESMF_DATEMOD 0     c   a   ESMF_DATE%CALENDAR+ESMF_DATEMOD /   ~        ESMF_CALENDAR+ESMF_CALENDARMOD 9   �  H   %   ESMF_CALENDAR%TYPE+ESMF_CALENDARMOD=TYPE 7   E  �   %   ESMF_CALENDAR%DIM+ESMF_CALENDARMOD=DIM K   �  �   %   ESMF_CALENDAR%DIMRUNNINGSUM+ESMF_CALENDARMOD=DIMRUNNINGSUM 7   }	  H   %   ESMF_CALENDAR%DIY+ESMF_CALENDARMOD=DIY 1   �	  H   %   ESMF_DATE%YEAR+ESMF_DATEMOD=YEAR 3   
  H   %   ESMF_DATE%MONTH+ESMF_DATEMOD=MONTH /   U
  H   %   ESMF_DATE%DAY+ESMF_DATEMOD=DAY +   �
  ^   a   ESMF_DATE%TOD+ESMF_DATEMOD ;   �
  H   %   ESMF_DATE%JULIANDAY+ESMF_DATEMOD=JULIANDAY ;   C  H   %   ESMF_DATE%DAYOFYEAR+ESMF_DATEMOD=DAYOFYEAR -   �  �       ESMF_TIMEMGR+ESMF_TIMEMGRMOD 9   ;  H   %   ESMF_TIMEMGR%NSTEP+ESMF_TIMEMGRMOD=NSTEP 6   �  _   a   ESMF_TIMEMGR%STEPSIZE+ESMF_TIMEMGRMOD 7   �  _   a   ESMF_TIMEMGR%STARTDATE+ESMF_TIMEMGRMOD 6   A  _   a   ESMF_TIMEMGR%STOPDATE+ESMF_TIMEMGRMOD 6   �  _   a   ESMF_TIMEMGR%BASEDATE+ESMF_TIMEMGRMOD 6   �  _   a   ESMF_TIMEMGR%CURRDATE+ESMF_TIMEMGRMOD 6   ^  _   a   ESMF_TIMEMGR%PREVDATE+ESMF_TIMEMGRMOD $   �  @       BNDTVVOLC+FILENAMES "   �  e       ENDRUN+ABORTUTILS 2   b  @      ENDRUN%PRESENT+ABORTUTILS=PRESENT &   �  L   e   ENDRUN%MSG+ABORTUTILS (   �  |       VALIDFACTORS+TIMEINTERP 0   j  <      VALIDFACTORS%ABS+TIMEINTERP=ABS .   �  @   e   VALIDFACTORS%FACT1+TIMEINTERP .   �  @   e   VALIDFACTORS%FACT2+TIMEINTERP #   &  �       READ_VOLCANIC_MASS P           READ_VOLCANIC_MASS%ESMF_CALENDAR+ESMF_CALENDARMOD=ESMF_CALENDAR L   �  H   %   READ_VOLCANIC_MASS%ESMF_CALENDAR%TYPE+ESMF_CALENDARMOD=TYPE J   �  �   %   READ_VOLCANIC_MASS%ESMF_CALENDAR%DIM+ESMF_CALENDARMOD=DIM ^   �  �   %   READ_VOLCANIC_MASS%ESMF_CALENDAR%DIMRUNNINGSUM+ESMF_CALENDARMOD=DIMRUNNINGSUM J     H   %   READ_VOLCANIC_MASS%ESMF_CALENDAR%DIY+ESMF_CALENDARMOD=DIY D   e  �      READ_VOLCANIC_MASS%ESMF_DATE+ESMF_DATEMOD=ESMF_DATE C     v   a   READ_VOLCANIC_MASS%ESMF_DATE%CALENDAR+ESMF_DATEMOD D   ~  H   %   READ_VOLCANIC_MASS%ESMF_DATE%YEAR+ESMF_DATEMOD=YEAR F   �  H   %   READ_VOLCANIC_MASS%ESMF_DATE%MONTH+ESMF_DATEMOD=MONTH B     H   %   READ_VOLCANIC_MASS%ESMF_DATE%DAY+ESMF_DATEMOD=DAY >   V  q   a   READ_VOLCANIC_MASS%ESMF_DATE%TOD+ESMF_DATEMOD A   �  m      READ_VOLCANIC_MASS%ESMF_TOD+ESMF_TODMOD=ESMF_TOD B   4  H   %   READ_VOLCANIC_MASS%ESMF_TOD%TYPE+ESMF_TODMOD=TYPE @   |  H   %   READ_VOLCANIC_MASS%ESMF_TOD%SEC+ESMF_TODMOD=SEC B   �  H   %   READ_VOLCANIC_MASS%ESMF_TOD%MSEC+ESMF_TODMOD=MSEC N     H   %   READ_VOLCANIC_MASS%ESMF_DATE%JULIANDAY+ESMF_DATEMOD=JULIANDAY N   T  H   %   READ_VOLCANIC_MASS%ESMF_DATE%DAYOFYEAR+ESMF_DATEMOD=DAYOFYEAR D   �  b      READ_VOLCANIC_MASS%ESMF_TIME+ESMF_TIMEMOD=ESMF_TIME B   �  H   %   READ_VOLCANIC_MASS%ESMF_TIME%DAY+ESMF_TIMEMOD=DAY >   F  q   a   READ_VOLCANIC_MASS%ESMF_TIME%TOD+ESMF_TIMEMOD M   �  �      READ_VOLCANIC_MASS%ESMF_TIMEMGR+ESMF_TIMEMGRMOD=ESMF_TIMEMGR L   g  H   %   READ_VOLCANIC_MASS%ESMF_TIMEMGR%NSTEP+ESMF_TIMEMGRMOD=NSTEP I   �  r   a   READ_VOLCANIC_MASS%ESMF_TIMEMGR%STEPSIZE+ESMF_TIMEMGRMOD J   !  r   a   READ_VOLCANIC_MASS%ESMF_TIMEMGR%STARTDATE+ESMF_TIMEMGRMOD I   �  r   a   READ_VOLCANIC_MASS%ESMF_TIMEMGR%STOPDATE+ESMF_TIMEMGRMOD I     r   a   READ_VOLCANIC_MASS%ESMF_TIMEMGR%BASEDATE+ESMF_TIMEMGRMOD I   w  r   a   READ_VOLCANIC_MASS%ESMF_TIMEMGR%CURRDATE+ESMF_TIMEMGRMOD I   �  r   a   READ_VOLCANIC_MASS%ESMF_TIMEMGR%PREVDATE+ESMF_TIMEMGRMOD "   [        GET_VOLCANIC_MASS O   f        GET_VOLCANIC_MASS%ESMF_CALENDAR+ESMF_CALENDARMOD=ESMF_CALENDAR K   �  H   %   GET_VOLCANIC_MASS%ESMF_CALENDAR%TYPE+ESMF_CALENDARMOD=TYPE I   -  �   %   GET_VOLCANIC_MASS%ESMF_CALENDAR%DIM+ESMF_CALENDARMOD=DIM ]   �  �   %   GET_VOLCANIC_MASS%ESMF_CALENDAR%DIMRUNNINGSUM+ESMF_CALENDARMOD=DIMRUNNINGSUM I   e   H   %   GET_VOLCANIC_MASS%ESMF_CALENDAR%DIY+ESMF_CALENDARMOD=DIY C   �   �      GET_VOLCANIC_MASS%ESMF_DATE+ESMF_DATEMOD=ESMF_DATE B   P!  u   a   GET_VOLCANIC_MASS%ESMF_DATE%CALENDAR+ESMF_DATEMOD C   �!  H   %   GET_VOLCANIC_MASS%ESMF_DATE%YEAR+ESMF_DATEMOD=YEAR E   "  H   %   GET_VOLCANIC_MASS%ESMF_DATE%MONTH+ESMF_DATEMOD=MONTH A   U"  H   %   GET_VOLCANIC_MASS%ESMF_DATE%DAY+ESMF_DATEMOD=DAY =   �"  p   a   GET_VOLCANIC_MASS%ESMF_DATE%TOD+ESMF_DATEMOD @   #  m      GET_VOLCANIC_MASS%ESMF_TOD+ESMF_TODMOD=ESMF_TOD A   z#  H   %   GET_VOLCANIC_MASS%ESMF_TOD%TYPE+ESMF_TODMOD=TYPE ?   �#  H   %   GET_VOLCANIC_MASS%ESMF_TOD%SEC+ESMF_TODMOD=SEC A   
$  H   %   GET_VOLCANIC_MASS%ESMF_TOD%MSEC+ESMF_TODMOD=MSEC M   R$  H   %   GET_VOLCANIC_MASS%ESMF_DATE%JULIANDAY+ESMF_DATEMOD=JULIANDAY M   �$  H   %   GET_VOLCANIC_MASS%ESMF_DATE%DAYOFYEAR+ESMF_DATEMOD=DAYOFYEAR C   �$  b      GET_VOLCANIC_MASS%ESMF_TIME+ESMF_TIMEMOD=ESMF_TIME A   D%  H   %   GET_VOLCANIC_MASS%ESMF_TIME%DAY+ESMF_TIMEMOD=DAY =   �%  p   a   GET_VOLCANIC_MASS%ESMF_TIME%TOD+ESMF_TIMEMOD L   �%  �      GET_VOLCANIC_MASS%ESMF_TIMEMGR+ESMF_TIMEMGRMOD=ESMF_TIMEMGR K   �&  H   %   GET_VOLCANIC_MASS%ESMF_TIMEMGR%NSTEP+ESMF_TIMEMGRMOD=NSTEP H   �&  q   a   GET_VOLCANIC_MASS%ESMF_TIMEMGR%STEPSIZE+ESMF_TIMEMGRMOD I   e'  q   a   GET_VOLCANIC_MASS%ESMF_TIMEMGR%STARTDATE+ESMF_TIMEMGRMOD H   �'  q   a   GET_VOLCANIC_MASS%ESMF_TIMEMGR%STOPDATE+ESMF_TIMEMGRMOD H   G(  q   a   GET_VOLCANIC_MASS%ESMF_TIMEMGR%BASEDATE+ESMF_TIMEMGRMOD H   �(  q   a   GET_VOLCANIC_MASS%ESMF_TIMEMGR%CURRDATE+ESMF_TIMEMGRMOD H   ))  q   a   GET_VOLCANIC_MASS%ESMF_TIMEMGR%PREVDATE+ESMF_TIMEMGRMOD (   �)  @   a   GET_VOLCANIC_MASS%LCHNK *   �)  �   a   GET_VOLCANIC_MASS%AEROSOL $   ~*  �       VOLCANIC_INITIALIZE Q   {+        VOLCANIC_INITIALIZE%ESMF_CALENDAR+ESMF_CALENDARMOD=ESMF_CALENDAR M   �+  H   %   VOLCANIC_INITIALIZE%ESMF_CALENDAR%TYPE+ESMF_CALENDARMOD=TYPE K   B,  �   %   VOLCANIC_INITIALIZE%ESMF_CALENDAR%DIM+ESMF_CALENDARMOD=DIM _   �,  �   %   VOLCANIC_INITIALIZE%ESMF_CALENDAR%DIMRUNNINGSUM+ESMF_CALENDARMOD=DIMRUNNINGSUM K   z-  H   %   VOLCANIC_INITIALIZE%ESMF_CALENDAR%DIY+ESMF_CALENDARMOD=DIY E   �-  �      VOLCANIC_INITIALIZE%ESMF_DATE+ESMF_DATEMOD=ESMF_DATE D   e.  w   a   VOLCANIC_INITIALIZE%ESMF_DATE%CALENDAR+ESMF_DATEMOD E   �.  H   %   VOLCANIC_INITIALIZE%ESMF_DATE%YEAR+ESMF_DATEMOD=YEAR G   $/  H   %   VOLCANIC_INITIALIZE%ESMF_DATE%MONTH+ESMF_DATEMOD=MONTH C   l/  H   %   VOLCANIC_INITIALIZE%ESMF_DATE%DAY+ESMF_DATEMOD=DAY ?   �/  r   a   VOLCANIC_INITIALIZE%ESMF_DATE%TOD+ESMF_DATEMOD B   &0  m      VOLCANIC_INITIALIZE%ESMF_TOD+ESMF_TODMOD=ESMF_TOD C   �0  H   %   VOLCANIC_INITIALIZE%ESMF_TOD%TYPE+ESMF_TODMOD=TYPE A   �0  H   %   VOLCANIC_INITIALIZE%ESMF_TOD%SEC+ESMF_TODMOD=SEC C   #1  H   %   VOLCANIC_INITIALIZE%ESMF_TOD%MSEC+ESMF_TODMOD=MSEC O   k1  H   %   VOLCANIC_INITIALIZE%ESMF_DATE%JULIANDAY+ESMF_DATEMOD=JULIANDAY O   �1  H   %   VOLCANIC_INITIALIZE%ESMF_DATE%DAYOFYEAR+ESMF_DATEMOD=DAYOFYEAR E   �1  b      VOLCANIC_INITIALIZE%ESMF_TIME+ESMF_TIMEMOD=ESMF_TIME C   ]2  H   %   VOLCANIC_INITIALIZE%ESMF_TIME%DAY+ESMF_TIMEMOD=DAY ?   �2  r   a   VOLCANIC_INITIALIZE%ESMF_TIME%TOD+ESMF_TIMEMOD N   3  �      VOLCANIC_INITIALIZE%ESMF_TIMEMGR+ESMF_TIMEMGRMOD=ESMF_TIMEMGR M   �3  H   %   VOLCANIC_INITIALIZE%ESMF_TIMEMGR%NSTEP+ESMF_TIMEMGRMOD=NSTEP J   4  s   a   VOLCANIC_INITIALIZE%ESMF_TIMEMGR%STEPSIZE+ESMF_TIMEMGRMOD K   �4  s   a   VOLCANIC_INITIALIZE%ESMF_TIMEMGR%STARTDATE+ESMF_TIMEMGRMOD J   �4  s   a   VOLCANIC_INITIALIZE%ESMF_TIMEMGR%STOPDATE+ESMF_TIMEMGRMOD J   h5  s   a   VOLCANIC_INITIALIZE%ESMF_TIMEMGR%BASEDATE+ESMF_TIMEMGRMOD J   �5  s   a   VOLCANIC_INITIALIZE%ESMF_TIMEMGR%CURRDATE+ESMF_TIMEMGRMOD J   N6  s   a   VOLCANIC_INITIALIZE%ESMF_TIMEMGR%PREVDATE+ESMF_TIMEMGRMOD 
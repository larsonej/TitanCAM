	  H  7   k820309              12.0        ÓlBO                                                                                                           
       /home/larsonej/titancam/src/models/lnd/clm2/src/biogeophys/SurfaceAlbedoMod.F90 SURFACEALBEDOMOD              SNOWALBEDO SOILALBEDO TWOSTREAM                      @                              
       SHR_SYS_FLUSH #         @                                                     #UNIT                @                                         #         @                                                    #SURFACEALBEDO%EXP    #SURFACEALBEDO%SQRT    #SURFACEALBEDO%MAX    #LBG    #UBG 	   #LBC 
   #UBC    #LBP    #UBP    #FRAC_DAY    #DAY_IN_YEAR    #ECCEN    #OBLIQR    #LAMBM0    #MVELPP                                                                                                    EXP                                                 SQRT                                                 MAX           
                                                       
                                  	                     
  @                               
                     
  @                                                    
  @                                                    
  @                                                    
  @                                   
                
  @                                   
                
  @                                   
                
  @                                   
                
  @                                   
                
  @                                   
      #         @                                                   #SNOWAGE%ESMF_CALENDAR    #SNOWAGE%ESMF_DATE    #SNOWAGE%ESMF_TOD     #SNOWAGE%ESMF_TIME &   #SNOWAGE%ESMF_TIMEMGR )   #SNOWAGE%MIN 1   #SNOWAGE%EXP 2   #SNOWAGE%MAX 3   #LBC 4   #UBC 5                                                       @                              'à                    #TYPE    #DIM    #DIMRUNNINGSUM    #DIY                  D                                                             D                                                           p          p            p                                        D                                        p                   p          p            p                                        D                                 Ø                            @                              '                    #CALENDAR    #YEAR    #MONTH    #DAY    #TOD    #JULIANDAY $   #DAYOFYEAR %                                                   à                      #SNOWAGE%ESMF_CALENDAR                 D                                 à                          D                                 è                          D                                 ð                                                                    ø              #SNOWAGE%ESMF_TOD                     @                               '                    #TYPE !   #SEC "   #MSEC #                 D                            !                                 D                            "                                D                            #                               D                            $                              D                            %                                @                         &     '                     #DAY '   #TOD (                D                            '                                                              (                          #SNOWAGE%ESMF_TOD                    @                          )     'È                   #NSTEP *   #STEPSIZE +   #STARTDATE ,   #STOPDATE -   #BASEDATE .   #CURRDATE /   #PREVDATE 0                D                            *                                                              +                           #SNOWAGE%ESMF_TIME &                                              ,            (              #SNOWAGE%ESMF_DATE                                               -            H             #SNOWAGE%ESMF_DATE                                               .            h             #SNOWAGE%ESMF_DATE                                               /                         #SNOWAGE%ESMF_DATE                                               0            ¨             #SNOWAGE%ESMF_DATE                                               1     MIN                                            2     EXP                                            3     MAX           
                                  4                     
                                  5                  i      fn#fn &   	  0   b   uapp(SURFACEALBEDOMOD    9  N   J  SHR_SYS_MOD *     R       SHR_SYS_FLUSH+SHR_SYS_MOD /   Ù  @   e   SHR_SYS_FLUSH%UNIT+SHR_SYS_MOD      B      SURFACEALBEDO "   [  <      SURFACEALBEDO%EXP #     =      SURFACEALBEDO%SQRT "   Ô  <      SURFACEALBEDO%MAX "     @   a   SURFACEALBEDO%LBG "   P  @   a   SURFACEALBEDO%UBG "     @   a   SURFACEALBEDO%LBC "   Ð  @   a   SURFACEALBEDO%UBC "     @   a   SURFACEALBEDO%LBP "   P  @   a   SURFACEALBEDO%UBP '     @   a   SURFACEALBEDO%FRAC_DAY *   Ð  @   a   SURFACEALBEDO%DAY_IN_YEAR $     @   a   SURFACEALBEDO%ECCEN %   P  @   a   SURFACEALBEDO%OBLIQR %     @   a   SURFACEALBEDO%LAMBM0 %   Ð  @   a   SURFACEALBEDO%MVELPP      *      SNOWAGE E   :        SNOWAGE%ESMF_CALENDAR+ESMF_CALENDARMOD=ESMF_CALENDAR A   ¹  H   %   SNOWAGE%ESMF_CALENDAR%TYPE+ESMF_CALENDARMOD=TYPE ?   	     %   SNOWAGE%ESMF_CALENDAR%DIM+ESMF_CALENDARMOD=DIM S   	     %   SNOWAGE%ESMF_CALENDAR%DIMRUNNINGSUM+ESMF_CALENDARMOD=DIMRUNNINGSUM ?   9
  H   %   SNOWAGE%ESMF_CALENDAR%DIY+ESMF_CALENDARMOD=DIY 9   
  £      SNOWAGE%ESMF_DATE+ESMF_DATEMOD=ESMF_DATE 8   $  k   a   SNOWAGE%ESMF_DATE%CALENDAR+ESMF_DATEMOD 9     H   %   SNOWAGE%ESMF_DATE%YEAR+ESMF_DATEMOD=YEAR ;   ×  H   %   SNOWAGE%ESMF_DATE%MONTH+ESMF_DATEMOD=MONTH 7     H   %   SNOWAGE%ESMF_DATE%DAY+ESMF_DATEMOD=DAY 3   g  f   a   SNOWAGE%ESMF_DATE%TOD+ESMF_DATEMOD 6   Í  m      SNOWAGE%ESMF_TOD+ESMF_TODMOD=ESMF_TOD 7   :  H   %   SNOWAGE%ESMF_TOD%TYPE+ESMF_TODMOD=TYPE 5     H   %   SNOWAGE%ESMF_TOD%SEC+ESMF_TODMOD=SEC 7   Ê  H   %   SNOWAGE%ESMF_TOD%MSEC+ESMF_TODMOD=MSEC C     H   %   SNOWAGE%ESMF_DATE%JULIANDAY+ESMF_DATEMOD=JULIANDAY C   Z  H   %   SNOWAGE%ESMF_DATE%DAYOFYEAR+ESMF_DATEMOD=DAYOFYEAR 9   ¢  b      SNOWAGE%ESMF_TIME+ESMF_TIMEMOD=ESMF_TIME 7     H   %   SNOWAGE%ESMF_TIME%DAY+ESMF_TIMEMOD=DAY 3   L  f   a   SNOWAGE%ESMF_TIME%TOD+ESMF_TIMEMOD B   ²  °      SNOWAGE%ESMF_TIMEMGR+ESMF_TIMEMGRMOD=ESMF_TIMEMGR A   b  H   %   SNOWAGE%ESMF_TIMEMGR%NSTEP+ESMF_TIMEMGRMOD=NSTEP >   ª  g   a   SNOWAGE%ESMF_TIMEMGR%STEPSIZE+ESMF_TIMEMGRMOD ?     g   a   SNOWAGE%ESMF_TIMEMGR%STARTDATE+ESMF_TIMEMGRMOD >   x  g   a   SNOWAGE%ESMF_TIMEMGR%STOPDATE+ESMF_TIMEMGRMOD >   ß  g   a   SNOWAGE%ESMF_TIMEMGR%BASEDATE+ESMF_TIMEMGRMOD >   F  g   a   SNOWAGE%ESMF_TIMEMGR%CURRDATE+ESMF_TIMEMGRMOD >   ­  g   a   SNOWAGE%ESMF_TIMEMGR%PREVDATE+ESMF_TIMEMGRMOD      <      SNOWAGE%MIN    P  <      SNOWAGE%EXP      <      SNOWAGE%MAX    È  @   a   SNOWAGE%LBC      @   a   SNOWAGE%UBC 
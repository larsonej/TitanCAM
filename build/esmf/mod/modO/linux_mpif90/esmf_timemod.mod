	  ¼  H   k820309              12.0        ªlBO                                                                                                           
       /home/larsonej/titancam/src/models/utils/esmf/src/Infrastructure/TimeMgmt/ESMF_TimeMod.F ESMF_TIMEMOD                                                    
                                                          
                                                              u #ESMF_TIMEINITIS    #ESMF_TIMEINITUNDEFINED    #ESMF_TIMECOPYINIT                                                           u #ESMF_TIMESETIS                                                           u #ESMF_TIMEGETIS                                                           u #ESMF_TIMEINCREMENTIS                                                           u #ESMF_TIMEDECREMENTIS 	                                                
                                                       0                @                                '                     #DAY    #TOD                  D                                                                                                                      #ESMF_TOD                    @                                '                    #TYPE    #SEC    #MSEC                  D                                                              D                                                             D                                               &         @@   X                                                        #ESMF_TIMEINITIS%PRESENT    #DAYS    #SECONDS    #RC    #ESMF_TIME                                                    PRESENT           
@ @                                                    
@ @                                                    F @                                           &         @@   X                                                        #ESMF_TIMEINITUNDEFINED%PRESENT    #RC    #ESMF_TIME                                                    PRESENT           F @                                           &         @@   X                                                        #ESMF_TIMECOPYINIT%PRESENT    #ORIG    #RC    #ESMF_TIME                                                    PRESENT           
@ @                                                   #ESMF_TIME              F @                                           #         @      X                                              #ESMF_TIMESETIS%PRESENT    #TIME    #DAYS    #SECONDS    #RC                                                    PRESENT           D @                                                    #ESMF_TIME              
@ @                                                    
@ @                                                    F @                                           #         @      X                                              #ESMF_TIMEGETIS%PRESENT     #TIME !   #DAYS "   #SECONDS #   #RC $                                                    PRESENT           
@ @                               !                    #ESMF_TIME              D @                               "                      D @                               #                      F @                               $            &         @@   X                                                        #ESMF_TIMEINCREMENTIS%PRESENT %   #TIME &   #DAYS '   #SECONDS (   #RC )   #ESMF_TIME                                               %     PRESENT           
@ @                               &                    #ESMF_TIME              
@ @                               '                     
@ @                               (                     F @                               )            &         @@   X                             	                           #ESMF_TIMEDECREMENTIS%PRESENT *   #TIME +   #DAYS ,   #SECONDS -   #RC .   #ESMF_TIME                                               *     PRESENT           
@ @                               +                    #ESMF_TIME              
@ @                               ,                     
@ @                               -                     F @                               .            %         @@                               /                   
       #ESMF_TIMEGETDAYS%PRESENT 0   #TIME 1   #RC 2                                              0     PRESENT           
@ @                               1                    #ESMF_TIME              F @                               2            #         @                                  3                  #ESMF_TIMECOPY%PRESENT 4   #TIME 5   #ORIG 6   #RC 7                                              4     PRESENT           D @                               5                     #ESMF_TIME              
@ @                               6                    #ESMF_TIME              F @                               7            #         @                                  8                  #ESMF_TIMEDIFF%PRESENT 9   #EARLYTIME :   #LATETIME ;   #DIFF <   #ISLATER =   #RC >                                              9     PRESENT           
@ @                               :                    #ESMF_TIME              
@ @                               ;                    #ESMF_TIME              D @                               <                     #ESMF_TIME              D @                               =                      F @                               >            #         @                                  ?                  #ESMF_TIMEPRINT%PRESENT @   #TIME A   #RC B                                              @     PRESENT           
@ @                               A                    #ESMF_TIME              F @                               B                   n      fn#fn "     @   J   ESMF_BASICUTILMOD    N  @   J   ESMF_TODMOD "            gen@ESMF_TIMEINIT !     T       gen@ESMF_TIMESET !   j  T       gen@ESMF_TIMEGET '   ¾  Z       gen@ESMF_TIMEINCREMENT '     Z       gen@ESMF_TIMEDECREMENT /   r  q       ESMF_SUCCESS+ESMF_BASICUTILMOD    ã  b       ESMF_TIME    E  H   !   ESMF_TIME%DAY      ^   a   ESMF_TIME%TOD %   ë  m       ESMF_TOD+ESMF_TODMOD *   X  H   %   ESMF_TOD%TYPE+ESMF_TODMOD )      H   %   ESMF_TOD%SEC+ESMF_TODMOD *   è  H   %   ESMF_TOD%MSEC+ESMF_TODMOD     0         ESMF_TIMEINITIS (   Ë  @      ESMF_TIMEINITIS%PRESENT %     @   a   ESMF_TIMEINITIS%DAYS (   K  @   a   ESMF_TIMEINITIS%SECONDS #     @   a   ESMF_TIMEINITIS%RC '   Ë         ESMF_TIMEINITUNDEFINED /   V  @      ESMF_TIMEINITUNDEFINED%PRESENT *     @   a   ESMF_TIMEINITUNDEFINED%RC "   Ö         ESMF_TIMECOPYINIT *   f	  @      ESMF_TIMECOPYINIT%PRESENT '   ¦	  W   a   ESMF_TIMECOPYINIT%ORIG %   ý	  @   a   ESMF_TIMECOPYINIT%RC    =
         ESMF_TIMESETIS '   Ê
  @      ESMF_TIMESETIS%PRESENT $   
  W   a   ESMF_TIMESETIS%TIME $   a  @   a   ESMF_TIMESETIS%DAYS '   ¡  @   a   ESMF_TIMESETIS%SECONDS "   á  @   a   ESMF_TIMESETIS%RC    !         ESMF_TIMEGETIS '   ®  @      ESMF_TIMEGETIS%PRESENT $   î  W   a   ESMF_TIMEGETIS%TIME $   E  @   a   ESMF_TIMEGETIS%DAYS '     @   a   ESMF_TIMEGETIS%SECONDS "   Å  @   a   ESMF_TIMEGETIS%RC %     ª       ESMF_TIMEINCREMENTIS -   ¯  @      ESMF_TIMEINCREMENTIS%PRESENT *   ï  W   a   ESMF_TIMEINCREMENTIS%TIME *   F  @   a   ESMF_TIMEINCREMENTIS%DAYS -     @   a   ESMF_TIMEINCREMENTIS%SECONDS (   Æ  @   a   ESMF_TIMEINCREMENTIS%RC %     ª       ESMF_TIMEDECREMENTIS -   °  @      ESMF_TIMEDECREMENTIS%PRESENT *   ð  W   a   ESMF_TIMEDECREMENTIS%TIME *   G  @   a   ESMF_TIMEDECREMENTIS%DAYS -     @   a   ESMF_TIMEDECREMENTIS%SECONDS (   Ç  @   a   ESMF_TIMEDECREMENTIS%RC !            ESMF_TIMEGETDAYS )     @      ESMF_TIMEGETDAYS%PRESENT &   Ç  W   a   ESMF_TIMEGETDAYS%TIME $     @   a   ESMF_TIMEGETDAYS%RC    ^         ESMF_TIMECOPY &   Ý  @      ESMF_TIMECOPY%PRESENT #     W   a   ESMF_TIMECOPY%TIME #   t  W   a   ESMF_TIMECOPY%ORIG !   Ë  @   a   ESMF_TIMECOPY%RC             ESMF_TIMEDIFF &   ª  @      ESMF_TIMEDIFF%PRESENT (   ê  W   a   ESMF_TIMEDIFF%EARLYTIME '   A  W   a   ESMF_TIMEDIFF%LATETIME #     W   a   ESMF_TIMEDIFF%DIFF &   ï  @   a   ESMF_TIMEDIFF%ISLATER !   /  @   a   ESMF_TIMEDIFF%RC    o  v       ESMF_TIMEPRINT '   å  @      ESMF_TIMEPRINT%PRESENT $   %  W   a   ESMF_TIMEPRINT%TIME "   |  @   a   ESMF_TIMEPRINT%RC 
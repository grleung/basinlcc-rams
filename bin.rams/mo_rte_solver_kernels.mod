  [_  `   k820309    w          19.1        +Òf                                                                                                          
       /swbuild/grleung1/bee-rams/src/radiate/rte-rrtmgp/rte/kernels/mo_rte_solver_kernels.F90 MO_RTE_SOLVER_KERNELS              LW_SOLVER_NOSCAT LW_SOLVER_2STREAM SW_SOLVER_NOSCAT SW_SOLVER_2STREAM                      @                              
                                                           
       WP WL                      @                              
       ZERO_ARRAY               @                                      u #ZERO_ARRAY_1D    #ZERO_ARRAY_2D    #ZERO_ARRAY_3D    #ZERO_ARRAY_4D    #         @     @                                               #NI    #ARRAY    n                      Ý              Czero_array_1D                                 
                                                                                                          	     p          5 O p            5 O p                          #         @     @                                              #NI    #NJ 	   #ARRAY 
   n                      á              Czero_array_2D                                 
                                                       
                                  	                                                    
                    	       p        5 O p        p          5 O p          5 O p            5 O p          5 O p                          #         @     @                                               #NI    #NJ    #NK    #ARRAY    n                      å              Czero_array_3D                                 
                                                       
                                                       
                                                                                                          	         p        5 O p        p        5 O p        p          5 O p          5 O p          5 O p            5 O p          5 O p          5 O p                          #         @     @                                               #NI    #NJ    #NK    #NL    #ARRAY    n                      é              Czero_array_4D                                 
                                                       
                                                       
                                                       
                                                                                                          	           p        5 O p        p        5 O p        p        5 O p        p          5 O p          5 O p          5 O p          5 O p            5 O p          5 O p          5 O p          5 O p                          #         @                                                      #NCOL    #NLAY    #NGPT    #TOP_AT_1    #NMUS    #DS    #WEIGHTS    #TAU    #LAY_SOURCE    #LEV_SOURCE_INC     #LEV_SOURCE_DEC !   #SFC_EMIS "   #SFC_SRC #   #INC_FLUX $   #FLUX_UP %   #FLUX_DN &   #DO_BROADBAND '   #BROADBAND_UP (   #BROADBAND_DN )   #DO_JACOBIANS *   #SFC_SRCJAC +   #FLUX_UPJAC ,   #DO_RESCALING -   #SSA .   #G /   n                                        Crte_lw_solver_noscat                               
  @                                                    
  @                                                    
  @                                                    
  @                                                    
                                                      
                                                     	        p        5  p        r    p        5  p        r    p          5  p        r      5  p        r      5  p        r        5  p        r      5  p        r      5  p        r                               
  @                                                  	     p          5  p        r        5  p        r                               
  @                                                  	 !       p        5  p        r    p        5  p        r    p          5  p        r      5  p        r      5  p        r        5  p        r      5  p        r      5  p        r                               
  @                                                  	 "       p        5  p        r    p        5  p        r    p          5  p        r      5  p        r      5  p        r        5  p        r      5  p        r      5  p        r                               
  @                                                   	 #       p        5  p        r    p        5  p        r    p          5  p        r      5  p        r      5  p        r        5  p        r      5  p        r      5  p        r                               
  @                              !                    	 $       p        5  p        r    p        5  p        r    p          5  p        r      5  p        r      5  p        r        5  p        r      5  p        r      5  p        r                               
  @                              "                    	 %     p        5  p        r    p          5  p        r      5  p        r        5  p        r      5  p        r                               
  @                              #                    	 &     p        5  p        r    p          5  p        r      5  p        r        5  p        r      5  p        r                               
  @                              $                    	 '     p        5  p        r    p          5  p        r      5  p        r        5  p        r      5  p        r                               D `                              %                    	 (        p         5  p        r    n                                           1p        5  p        r    p          5  p        r       5  p        r    n                                      1  5  p        r        5  p        r       5  p        r    n                                      1  5  p        r                                         D `                              &                    	 )        p         5  p        r    n                                           1p        5  p        r    p          5  p        r       5  p        r    n                                      1  5  p        r        5  p        r       5  p        r    n                                      1  5  p        r                                          
  @                               '                    D `                              (                    	 *      p        5  p        r    p          5  p        r       5  p        r    n                                       1    5  p        r       5  p        r    n                                      1                                    D `                              )                    	 +      p        5  p        r    p          5  p        r       5  p        r    n                                       1    5  p        r       5  p        r    n                                      1                                     
  @                               *                    
  @                              +                    	 ,     p        5  p        r    p          5  p        r      5  p        r        5  p        r      5  p        r                               D `                              ,                    	 -      p        5  p        r    p          5  p        r       5  p        r    n                                       1    5  p        r       5  p        r    n                                      1                                     
  @                               -                    
  @                              .                    	 .       p        5  p        r    p        5  p        r    p          5  p        r      5  p        r      5  p        r        5  p        r      5  p        r      5  p        r                               
  @                              /                    	 /       p        5  p        r    p        5  p        r    p          5  p        r      5  p        r      5  p        r        5  p        r      5  p        r      5  p        r                      #         @                                  0                    #NCOL 1   #NLAY 2   #NGPT 3   #TOP_AT_1 4   #TAU 5   #SSA 6   #G 7   #LAY_SOURCE 8   #LEV_SOURCE_INC 9   #LEV_SOURCE_DEC :   #SFC_EMIS ;   #SFC_SRC <   #INC_FLUX =   #FLUX_UP >   #FLUX_DN ?   n                                         Crte_lw_solver_2stream                              
  @                               1                     
  @                               2                     
                                  3                     
  @                               4                    
                                 5                    	 5       p        5  p        r 2   p        5  p        r 1   p          5  p        r 1     5  p        r 2     5  p        r 3       5  p        r 1     5  p        r 2     5  p        r 3                              
                                 6                    	 6       p        5  p        r 2   p        5  p        r 1   p          5  p        r 1     5  p        r 2     5  p        r 3       5  p        r 1     5  p        r 2     5  p        r 3                              
                                 7                    	 7       p        5  p        r 2   p        5  p        r 1   p          5  p        r 1     5  p        r 2     5  p        r 3       5  p        r 1     5  p        r 2     5  p        r 3                              
                                 8                    	 8       p        5  p        r 2   p        5  p        r 1   p          5  p        r 1     5  p        r 2     5  p        r 3       5  p        r 1     5  p        r 2     5  p        r 3                              
                                 9                    	 9       p        5  p        r 2   p        5  p        r 1   p          5  p        r 1     5  p        r 2     5  p        r 3       5  p        r 1     5  p        r 2     5  p        r 3                              
                                 :                    	 :       p        5  p        r 2   p        5  p        r 1   p          5  p        r 1     5  p        r 2     5  p        r 3       5  p        r 1     5  p        r 2     5  p        r 3                              
                                 ;                    	 ;     p        5  p        r 1   p          5  p        r 1     5  p        r 3       5  p        r 1     5  p        r 3                              
                                 <                    	 <     p        5  p        r 1   p          5  p        r 1     5  p        r 3       5  p        r 1     5  p        r 3                              
                                 =                    	 =     p        5  p        r 1   p          5  p        r 1     5  p        r 3       5  p        r 1     5  p        r 3                              D @                              >                    	 >        p         5  p        r 2   n                                           1p        5  p        r 1   p          5  p        r 1      5  p        r 2   n                                      1  5  p        r 3       5  p        r 1      5  p        r 2   n                                      1  5  p        r 3                                        D @                              ?                    	 ?        p         5  p        r 2   n                                           1p        5  p        r 1   p          5  p        r 1      5  p        r 2   n                                      1  5  p        r 3       5  p        r 1      5  p        r 2   n                                      1  5  p        r 3                               #         @                                 @                    #NCOL A   #NLAY B   #NGPT C   #TOP_AT_1 D   #TAU E   #MU0 F   #INC_FLUX_DIR G   #FLUX_DIR H   n                                           Crte_sw_solver_noscat                            
                                  A                     
                                  B                     
                                  C                     
                                  D                    
                                 E                    	 I       p        5  p        r B   p        5  p        r A   p          5  p        r A     5  p        r B     5  p        r C       5  p        r A     5  p        r B     5  p        r C                              
                                 F                    	 J     p        5  p        r A   p          5  p        r A     5  p        r B       5  p        r A     5  p        r B                              
                                 G                    	 K     p        5  p        r A   p          5  p        r A     5  p        r C       5  p        r A     5  p        r C                              D                                H                    	 L        p         5  p        r B   n                                           1p        5  p        r A   p          5  p        r A      5  p        r B   n                                      1  5  p        r C       5  p        r A      5  p        r B   n                                      1  5  p        r C                               #         @                                  I                    #NCOL J   #NLAY K   #NGPT L   #TOP_AT_1 M   #TAU N   #SSA O   #G P   #MU0 Q   #SFC_ALB_DIR R   #SFC_ALB_DIF S   #INC_FLUX_DIR T   #FLUX_UP U   #FLUX_DN V   #FLUX_DIR W   #HAS_DIF_BC X   #INC_FLUX_DIF Y   #DO_BROADBAND Z   #BROADBAND_UP [   #BROADBAND_DN \   #BROADBAND_DIR ]   n                                            Crte_sw_solver_2stream                           
  @                               J                     
  @                               K                     
                                  L                     
  @                               M                    
                                 N                    	 M       p        5  p        r K   p        5  p        r J   p          5  p        r J     5  p        r K     5  p        r L       5  p        r J     5  p        r K     5  p        r L                              
                                 O                    	 N       p        5  p        r K   p        5  p        r J   p          5  p        r J     5  p        r K     5  p        r L       5  p        r J     5  p        r K     5  p        r L                              
                                 P                    	 O       p        5  p        r K   p        5  p        r J   p          5  p        r J     5  p        r K     5  p        r L       5  p        r J     5  p        r K     5  p        r L                              
  @                              Q                    	 P     p        5  p        r J   p          5  p        r J     5  p        r K       5  p        r J     5  p        r K                              
                                 R                    	 Q     p        5  p        r J   p          5  p        r J     5  p        r L       5  p        r J     5  p        r L                              
                                 S                    	 R     p        5  p        r J   p          5  p        r J     5  p        r L       5  p        r J     5  p        r L                              
                                 T                    	 S     p        5  p        r J   p          5  p        r J     5  p        r L       5  p        r J     5  p        r L                                                              U                    	 T        p         5  p        r K   n                                           1p        5  p        r J   p          5  p        r J      5  p        r K   n                                      1  5  p        r L       5  p        r J      5  p        r K   n                                      1  5  p        r L                                                                        V                    	 U        p         5  p        r K   n                                           1p        5  p        r J   p          5  p        r J      5  p        r K   n                                      1  5  p        r L       5  p        r J      5  p        r K   n                                      1  5  p        r L                                                                        W                    	 V        p         5  p        r K   n                                           1p        5  p        r J   p          5  p        r J      5  p        r K   n                                      1  5  p        r L       5  p        r J      5  p        r K   n                                      1  5  p        r L                                         
                                  X                    
                                 Y                    	 W     p        5  p        r J   p          5  p        r J     5  p        r L       5  p        r J     5  p        r L                               
                                  Z                    D @                              [                    	 X      p        5  p        r J   p          5  p        r J      5  p        r K   n                                       1    5  p        r J      5  p        r K   n                                      1                                    D @                              \                    	 Y      p        5  p        r J   p          5  p        r J      5  p        r K   n                                       1    5  p        r J      5  p        r K   n                                      1                                    D @                              ]                    	 Z      p        5  p        r J   p          5  p        r J      5  p        r K   n                                       1    5  p        r J      5  p        r K   n                                      1                                  v      fn#fn +     V   b   uapp(MO_RTE_SOLVER_KERNELS    l  @   J  ISO_C_BINDING    ¬  F   J  MO_RTE_KIND "   ò  K   J  MO_RTE_UTIL_ARRAY 1   =         gen@ZERO_ARRAY+MO_RTE_UTIL_ARRAY 0   É  ­      ZERO_ARRAY_1D+MO_RTE_UTIL_ARRAY 3   v  @   a   ZERO_ARRAY_1D%NI+MO_RTE_UTIL_ARRAY 6   ¶  ¤   a   ZERO_ARRAY_1D%ARRAY+MO_RTE_UTIL_ARRAY 0   Z  µ      ZERO_ARRAY_2D+MO_RTE_UTIL_ARRAY 3     @   a   ZERO_ARRAY_2D%NI+MO_RTE_UTIL_ARRAY 3   O  @   a   ZERO_ARRAY_2D%NJ+MO_RTE_UTIL_ARRAY 6     ü   a   ZERO_ARRAY_2D%ARRAY+MO_RTE_UTIL_ARRAY 0     ½      ZERO_ARRAY_3D+MO_RTE_UTIL_ARRAY 3   H  @   a   ZERO_ARRAY_3D%NI+MO_RTE_UTIL_ARRAY 3     @   a   ZERO_ARRAY_3D%NJ+MO_RTE_UTIL_ARRAY 3   È  @   a   ZERO_ARRAY_3D%NK+MO_RTE_UTIL_ARRAY 6     T  a   ZERO_ARRAY_3D%ARRAY+MO_RTE_UTIL_ARRAY 0   \	  Å      ZERO_ARRAY_4D+MO_RTE_UTIL_ARRAY 3   !
  @   a   ZERO_ARRAY_4D%NI+MO_RTE_UTIL_ARRAY 3   a
  @   a   ZERO_ARRAY_4D%NJ+MO_RTE_UTIL_ARRAY 3   ¡
  @   a   ZERO_ARRAY_4D%NK+MO_RTE_UTIL_ARRAY 3   á
  @   a   ZERO_ARRAY_4D%NL+MO_RTE_UTIL_ARRAY 6   !  ¬  a   ZERO_ARRAY_4D%ARRAY+MO_RTE_UTIL_ARRAY !   Í  ú      LW_SOLVER_NOSCAT &   Ç  @   a   LW_SOLVER_NOSCAT%NCOL &     @   a   LW_SOLVER_NOSCAT%NLAY &   G  @   a   LW_SOLVER_NOSCAT%NGPT *     @   a   LW_SOLVER_NOSCAT%TOP_AT_1 &   Ç  @   a   LW_SOLVER_NOSCAT%NMUS $       a   LW_SOLVER_NOSCAT%DS )     ´   a   LW_SOLVER_NOSCAT%WEIGHTS %   O    a   LW_SOLVER_NOSCAT%TAU ,   ã    a   LW_SOLVER_NOSCAT%LAY_SOURCE 0   w    a   LW_SOLVER_NOSCAT%LEV_SOURCE_INC 0       a   LW_SOLVER_NOSCAT%LEV_SOURCE_DEC *     $  a   LW_SOLVER_NOSCAT%SFC_EMIS )   Ã  $  a   LW_SOLVER_NOSCAT%SFC_SRC *   ç  $  a   LW_SOLVER_NOSCAT%INC_FLUX )     ?  a   LW_SOLVER_NOSCAT%FLUX_UP )   J  ?  a   LW_SOLVER_NOSCAT%FLUX_DN .      @   a   LW_SOLVER_NOSCAT%DO_BROADBAND .   É     a   LW_SOLVER_NOSCAT%BROADBAND_UP .   _"    a   LW_SOLVER_NOSCAT%BROADBAND_DN .   õ#  @   a   LW_SOLVER_NOSCAT%DO_JACOBIANS ,   5$  $  a   LW_SOLVER_NOSCAT%SFC_SRCJAC ,   Y%    a   LW_SOLVER_NOSCAT%FLUX_UPJAC .   ï&  @   a   LW_SOLVER_NOSCAT%DO_RESCALING %   /'    a   LW_SOLVER_NOSCAT%SSA #   Ã(    a   LW_SOLVER_NOSCAT%G "   W*  b      LW_SOLVER_2STREAM '   ¹+  @   a   LW_SOLVER_2STREAM%NCOL '   ù+  @   a   LW_SOLVER_2STREAM%NLAY '   9,  @   a   LW_SOLVER_2STREAM%NGPT +   y,  @   a   LW_SOLVER_2STREAM%TOP_AT_1 &   ¹,    a   LW_SOLVER_2STREAM%TAU &   M.    a   LW_SOLVER_2STREAM%SSA $   á/    a   LW_SOLVER_2STREAM%G -   u1    a   LW_SOLVER_2STREAM%LAY_SOURCE 1   	3    a   LW_SOLVER_2STREAM%LEV_SOURCE_INC 1   4    a   LW_SOLVER_2STREAM%LEV_SOURCE_DEC +   16  $  a   LW_SOLVER_2STREAM%SFC_EMIS *   U7  $  a   LW_SOLVER_2STREAM%SFC_SRC +   y8  $  a   LW_SOLVER_2STREAM%INC_FLUX *   9  ?  a   LW_SOLVER_2STREAM%FLUX_UP *   Ü;  ?  a   LW_SOLVER_2STREAM%FLUX_DN !   >  ÿ       SW_SOLVER_NOSCAT &   ?  @   a   SW_SOLVER_NOSCAT%NCOL &   Z?  @   a   SW_SOLVER_NOSCAT%NLAY &   ?  @   a   SW_SOLVER_NOSCAT%NGPT *   Ú?  @   a   SW_SOLVER_NOSCAT%TOP_AT_1 %   @    a   SW_SOLVER_NOSCAT%TAU %   ®A  $  a   SW_SOLVER_NOSCAT%MU0 .   ÒB  $  a   SW_SOLVER_NOSCAT%INC_FLUX_DIR *   öC  ?  a   SW_SOLVER_NOSCAT%FLUX_DIR "   5F  ·      SW_SOLVER_2STREAM '   ìG  @   a   SW_SOLVER_2STREAM%NCOL '   ,H  @   a   SW_SOLVER_2STREAM%NLAY '   lH  @   a   SW_SOLVER_2STREAM%NGPT +   ¬H  @   a   SW_SOLVER_2STREAM%TOP_AT_1 &   ìH    a   SW_SOLVER_2STREAM%TAU &   J    a   SW_SOLVER_2STREAM%SSA $   L    a   SW_SOLVER_2STREAM%G &   ¨M  $  a   SW_SOLVER_2STREAM%MU0 .   ÌN  $  a   SW_SOLVER_2STREAM%SFC_ALB_DIR .   ðO  $  a   SW_SOLVER_2STREAM%SFC_ALB_DIF /   Q  $  a   SW_SOLVER_2STREAM%INC_FLUX_DIR *   8R  ?  a   SW_SOLVER_2STREAM%FLUX_UP *   wT  ?  a   SW_SOLVER_2STREAM%FLUX_DN +   ¶V  ?  a   SW_SOLVER_2STREAM%FLUX_DIR -   õX  @   a   SW_SOLVER_2STREAM%HAS_DIF_BC /   5Y  $  a   SW_SOLVER_2STREAM%INC_FLUX_DIF /   YZ  @   a   SW_SOLVER_2STREAM%DO_BROADBAND /   Z    a   SW_SOLVER_2STREAM%BROADBAND_UP /   /\    a   SW_SOLVER_2STREAM%BROADBAND_DN 0   Å]    a   SW_SOLVER_2STREAM%BROADBAND_DIR 
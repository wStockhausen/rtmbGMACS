## ———————————————————————————————————————————————————————————————————————————————————— ##
RECRUITMENT                                                                             ##
## ———————————————————————————————————————————————————————————————————————————————————— ##
   ival              lb        ub        phz    prior     p1      p2          # parameter
    8.0             -10        20         -2      0     -10.0    20.0         # logR0
    8.0    	        -10        20          1      0      10.0    20.0         # logRini, to estimate if NOT initialized at unfished (n68)
    8.0             -10        20          1      0      10.0    20.0         # logRbar, to estimate if NOT initialized at unfished
   -0.9    	        -10         0.75      -4      0     -10.0     0.75        # ln(sigma_R)
    0.75              0.20      1.00      -2      3       3.0     2.00        # steepness
    0.01    	        0.0001    1.00      -3      3       1.01    1.01        # recruitment autocorrelation
   32.5              10        42.5       -4      0      32.5     2.25        # recruitment size distribution expected value (males or combined)
    1.0               0.1      10	        -4      0       0.1     5.0         # recruitment sie distribution scale (variance component) (males or combined)
    0.0    	        -10        10         -4      0       0.0    20.00        # recruitment size distribution expected value (females)
    0.00   	        -10        10 	      -3      0       0.0    20.0         # recruitment size distribution scale (variance component) (females)
END TABLE
## rec devs
## sex ratio & devs
END SECTION #--RECRUITMENT
END CTL FILE

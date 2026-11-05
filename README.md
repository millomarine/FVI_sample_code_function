The file FVI_MC.R contains the function to calculate the FVI, according to Marín-Enríquez et al. The main arguments of the function are
K and Linf= the K and the L infinity parameters of the von Bertalanffy growth model, Tmax, the maximum observed age and L50, the size of first maturation at 50%.
The function also incorporates the uncertainity of each argument (Low and Up, for examle KLow and Kup for the K parameter), and the potential correlation between K 
and Linf (from -1 to 1) in the argument cor_K_Linf. The function issues a warning if the value of cor_K_Linf is positive (not a common case).

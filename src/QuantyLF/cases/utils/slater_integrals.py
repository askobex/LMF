def get_slater_integrals(ion, oxy):

    zeta_3d, F2dd, F4dd, zeta_2p, F2pd, G1pd, G3pd, Xzeta_3d, XF2dd, XF4dd = 0,0,0,0,0,0,0,0,0,0

    ## Slater integral values

    if ion == 29:        
        zeta_3d = 0.102
        F2dd = 12.854
        F4dd = 7.980 ## initial state parameters
        zeta_2p = 13.498
        F2pd = 8.177
        G1pd = 6.169
        G3pd = 3.510 ##  final  state parameters
        Xzeta_3d = 0.124
        XF2dd = 13.611
        XF4dd = 8.457 ##  final  state parameters

    elif ion == 28:
        if oxy == 2:
            zeta_3d = 0.083 ## Crispy
            F2dd = 12.234 ## Thesis de Groot
            F4dd = 7.598 ## Thesis de Groot
            zeta_2p = 11.507 ## Thesis de Groot
            F2pd = 7.721 ## Thesis de Groot
            G1pd = 5.787 ## Thesis de Groot
            G3pd = 3.291 ## Thesis de Groot
            Xzeta_3d = 0.102 ## Thesis de Groot
            XF2dd = 13.005
            XF4dd = 8.084
        elif oxy == 3:
            zeta_3d = 0.091 ## Crispy
            F2dd = 13.277 ## Thesis de Groot
            F4dd = 8.295 ## Thesis de Groot
            zeta_2p = 11.506 ## Thesis de Groot
            F2pd = 8.350 ## Thesis de Groot
            G1pd = 6.332 ## Thesis de Groot
            G3pd = 3.603 ## Thesis de Groot
            Xzeta_3d = 0.112 ## Thesis de Groot
            XF2dd = 14.022 ## Thesis de Groot
            XF4dd = 8.764 ## Thesis de Groot
        else:
            print("No data available for this ion and valence configuration...")
        

    elif ion == 27:
        if oxy == 2:
            zeta_3d = 0.066 ## Crispy
            F2dd = 11.605 ## Thesis de Groot
            F4dd = 7.209 ## Thesis de Groot
            zeta_2p = 9.746 ## Thesis de Groot
            F2pd = 7.260 ## Thesis de Groot
            G1pd = 5.397 ## Thesis de Groot
            G3pd = 3.069 ## Thesis de Groot
            Xzeta_3d = 0.092 ## Thesis de Groot
            XF2dd = 12.396 ## Thesis de Groot
            XF4dd = 7.708 ## Thesis de Groot
        elif oxy == 3:
            zeta_3d = 0.074 ## Crispy
            F2dd = 12.663 ## Thesis de Groot
            F4dd = 7.917 ## Thesis de Groot
            zeta_2p = 9.748 ## Thesis de Groot
            F2pd = 7.900 ## Thesis de Groot
            G1pd = 5.961 ## Thesis de Groot
            G3pd = 3.386 ## Thesis de Groot
            Xzeta_3d = 0.082 ## Thesis de Groot
            XF2dd = 13.422 ## Thesis de Groot
            XF4dd = 8.395 ## Thesis de Groot
        else:
            print("No data available for this ion and valence configuration...")
            
        

    elif ion == 26:
        if oxy == 2:
            zeta_3d = 0.052 ## Crispy
            F2dd = 10.966 ## Thesis de Groot
            F4dd = 6.815 ## Thesis de Groot
            zeta_2p = 8.200 ## Thesis de Groot
            F2pd = 6.793 ## Thesis de Groot
            G1pd = 5.004 ## Thesis de Groot
            G3pd = 2.844 ## Thesis de Groot
            Xzeta_3d = 0.067 ## Thesis de Groot
            XF2dd = 11.779 ## Thesis de Groot
            XF4dd = 7.327 ## Thesis de Groot
        elif oxy == 3:
            zeta_3d = 0.059 ## Crispy
            F2dd = 12.043 ## Thesis de Groot
            F4dd = 7.535 ## Thesis de Groot
            zeta_2p = 8.199 ## Thesis de Groot
            F2pd = 7.446 ## Thesis de Groot
            G1pd = 5.566 ## Thesis de Groot
            G3pd = 3.166 ## Thesis de Groot
            Xzeta_3d = 0.074 ## Thesis de Groot
            XF2dd = 12.818 ## Thesis de Groot
            XF4dd = 8.023 ## Thesis de Groot
        else:
            print("No data available for this ion and valence configuration...")
            
        

    elif ion == 25:
        if oxy == 2:
            zeta_3d = 0.040 ## Crispy
            F2dd = 10.316 ## Thesis de Groot
            F4dd = 6.414 ## Thesis de Groot
            zeta_2p = 6.846 ## Thesis de Groot
            F2pd = 6.321 ## Thesis de Groot
            G1pd = 4.606 ## Thesis de Groot
            G3pd = 2.618 ## Thesis de Groot
            Xzeta_3d = 0.053 ## Thesis de Groot
            XF2dd = 11.155 ## Thesis de Groot
            XF4dd = 6.943 ## Thesis de Groot
        elif oxy == 3:
            zeta_3d = 0.046 ## Crispy
            F2dd = 11.415 ## Thesis de Groot
            F4dd = 7.148 ## Thesis de Groot
            zeta_2p = 6.845 ## Thesis de Groot
            F2pd = 6.988 ## Thesis de Groot
            G1pd = 5.179 ## Thesis de Groot
            G3pd = 2.945 ## Thesis de Groot
            Xzeta_3d = 0.059 ## Thesis de Groot
            XF2dd = 12.210 ## Thesis de Groot
            XF4dd = 7.649 ## Thesis de Groot
        elif oxy == 4:
            zeta_3d = 0.052 ## Crispy
            F2dd = 12.416 ## Thesis de Groot
            F4dd = 7.820 ## Thesis de Groot
            zeta_2p = 6.845 ## Thesis de Groot
            F2pd = 7.658 ## Thesis de Groot
            G1pd = 5.776 ## Thesis de Groot
            G3pd = 3.288 ## Thesis de Groot
            Xzeta_3d = 0.066 ## Thesis de Groot
            XF2dd = 13.177 ## Thesis de Groot
            XF4dd = 8.299 ## Thesis de Groot
        else:
            print("No data available for this ion and valence configuration...")
            
        

    elif ion == 24:
        if oxy == 2:
            zeta_3d = 0.030 ## Crispy
            F2dd = 9.649 ## Thesis de Groot
            F4dd = 6.002 ## Thesis de Groot
            zeta_2p = 5.668 ## Thesis de Groot
            F2pd = 5.841 ## Thesis de Groot
            G1pd = 4.204 ## Thesis de Groot
            G3pd = 2.388 ## Thesis de Groot
            Xzeta_3d = 0.041 ## Thesis de Groot
            XF2dd = 10.522 ## Thesis de Groot
            XF4dd = 6.552 ## Thesis de Groot
        elif oxy == 3:
            zeta_3d = 0.035 ## Crispy
            F2dd = 10.777 ## Thesis de Groot
            F4dd = 6.755 ## Thesis de Groot
            zeta_2p = 5.667 ## Thesis de Groot
            F2pd = 6.526 ## Thesis de Groot
            G1pd = 4.788 ## Thesis de Groot
            G3pd = 2.722 ## Thesis de Groot
            Xzeta_3d = 0.047 ## Thesis de Groot
            XF2dd = 11.596 ## Thesis de Groot
            XF4dd = 7.270 ## Thesis de Groot
        else:
            print("No data available for this ion and valence configuration...")
            
        

     elif ion == 23: ## Crispy
        if oxy == 2:
            nd = 3
            zeta_3d = 0.022
            F2dd = 8.962
            F4dd = 5.577
            zeta_2p = 4.650
            F2pd = 5.351
            G1pd = 3.793
            G3pd  2.154
            Xzeta_3d = 0.031
            XF2dd = 9.876
            XF4dd = 6.153
        elif oxy == 3:   ## Crispy
            nd = 2
            zeta_3d = 0.027
            F2dd = 10.127
            F4dd = 6.354
            zeta_2p = 4.650
            F2pd = 6.056
            G1pd = 4.390
            G3pd = .496
            Xzeta_3d = 0.036
            XF2dd = 10.974
            XF4dd = 6.888
        elif oxy == 4:  ## Crispy
            nd = 1
            zeta_3d = 0.031
            F2dd = 0.000
            F4dd = 0.000
            zeta_2p = 4.650
            F2pd = 6.758
            G1pd = 5.012
            G3pd = 2.852
            Xzeta_3d = 0.042
            XF2dd = 11.965
            XF4dd = 7.554
        elif oxy == 5:  ## Crispy
            nd = 0
            zeta_3d = 0.000
            F2dd = 0.000
            F4dd = 0.000
            zeta_2p = 4.652
            F2pd = 7.460
            G1pd = 5.661
            G3pd = 3.226
            Xzeta_3d = 0.047
            XF2dd = 0.000
            XF4dd = 0.000
        else:
            print("No data available for this ion and valence configuration...")


    elif ion == 22:
        if oxy == 2:
            nd = 2
            zeta_3d=0.016
            F2dd=8.243
            F4dd=5.132
            zeta_2p=3.776
            F2pd=4.849
            G1pd=3.376
            G3pd=1.917
            Xzeta_3d=0.023
            XF2dd=9.213
            XF4dd=5.744
        elif oxy == 3:
            nd = 1
            zeta_3d=0.019
            F2dd=5.907
            F4dd=3.711
            zeta_2p=3.710
            F2pd=3.187
            G1pd=3.393
            G3pd=1.258
            Xzeta_3d=0.027
            XF2dd=10.342
            XF4dd=6.499
        elif oxy == 4:
            nd = 0
            zeta_3d=0.019
            F2dd=0.000
            F4dd=0.000
            zeta_2p=3.776
            F2pd=6.301
            G1pd=4.626
            G3pd=2.632
            Xzeta_3d=0.032
            XF2dd=0.000
            XF4dd=0.000
        else:
            print("No data available for this ion and valence configuration...")

    elif ion == 21:
        
        zeta_3d = 0.010
        F2dd = 0
        F4dd = 0
        zeta_2p = 3.032
        F2pd = 4.332
        G1pd = 2.950
        G3pd = 1.674
        Xzeta_3d = 0.017
        XF2dd = 8.530
        XF4dd = 5.321

    else:
        print("Could not recognize the ion name...")
        
    
    return zeta_3d, F2dd, F4dd, zeta_2p, F2pd, G1pd, G3pd, Xzeta_3d, XF2dd, XF4dd


import numpy as np


def dy_dt(y, t, IL2, IL15, k1fwd, k4fwd, k5rev, k6rev, k10rev, k11rev):
    # IL2 in nM
    IL2Ra = y[0]
    IL2Rb = y[1]
    gc = y[2]
    IL2_IL2Ra = y[3]
    IL2_IL2Rb = y[4]
    IL2_gc = y[5]
    IL2_IL2Ra_IL2Rb = y[6]
    IL2_IL2Ra_gc = y[7]
    IL2_IL2Rb_gc = y[8]
    IL2_IL2Ra_IL2Rb_gc = y[9]
    
    ###Ali's Addition
    IL15Ra = y[10]
    IL15_IL15Ra = y[11]
    IL15_IL2Rb = y[12]
    IL15_gc = y[13]
    IL15_IL15Ra_IL2Rb = y[14]
    IL15_IL15Ra_gc = y[15]
    IL15_IL2Rb_gc = y[16]
    IL15_IL15Ra_IL2Rb_gc = y[17]
    


    # The receptor-receptor forward rate is largely going to be determined by plasma membrane diffusion
    # so we're going to assume it's shared.
    k5fwd = k6fwd = k7fwd = k8fwd = k9fwd = k10fwd = k11fwd = k12fwd = k4fwd
    
    
    ###Ali's Addition
    k17fwd = k18fwd = k19fwd = k20fwd = k21fwd = k22fwd = k23fwd = k24fwd = k16fwd

    

    # These are probably measured in the literature
    k1rev = k1fwd * 10 # doi:10.1016/j.jmb.2004.04.038, 10 nM
    k2fwd = k1fwd
    k2rev = k2fwd * 144 # doi:10.1016/j.jmb.2004.04.038, 144 nM
    k3fwd = k1fwd / 10.0 # Very weak, > 50 uM. Voss, et al (1993). PNAS. 90, 2428–2432.
    k3rev = 50000 * k3fwd
    
    
    

    ###Ali's Addition to literature values for k values for IL-15





    # To satisfy detailed balance these relationships should hold
    # _Based on initial assembly steps
    k4rev = k1fwd * k4fwd * k6rev * k3rev / k1rev / k6fwd / k3fwd
    k7rev = k3fwd * k7fwd * k2rev * k5rev / k2fwd / k5fwd / k3rev
    k12rev = k2fwd * k12fwd * k1rev * k11rev / k11fwd / k1fwd / k2rev
    # _Based on formation of full complex
    k9rev = k2rev * k10rev * k12rev / k2fwd / k10fwd / k12fwd / k3rev / k6rev * k3fwd * k6fwd * k9fwd
    k8rev = k2rev * k10rev * k12rev / k2fwd / k10fwd / k12fwd / k7rev / k3rev * k3fwd * k7fwd * k8fwd



    ###Ali's Addition
    # To satisfy detailed balance these relationships should hold
    # _Based on initial assembly steps
    k16rev = k13fwd * k16fwd * k18rev * k15rev / k13rev / k18fwd / k15fwd
    k19rev = k15fwd * k19fwd * k14rev * k17rev / k14fwd / k17fwd / k15rev
    k24rev = k14fwd * k24fwd * k13rev * k23rev / k23fwd / k13fwd / k14rev
    # _Based on formation of full complex
    k21rev = k14rev * k22rev * k24rev / k14fwd / k22fwd / k24fwd / k15rev / k18rev * k15fwd * k18fwd * k21fwd
    k20rev = k14rev * k22rev * k24rev / k14fwd / k22fwd / k24fwd / k19rev / k15rev * k15fwd * k19fwd * k20fwd
    
    
    
    
    

    dydt = np.zeros(y.shape, dtype = np.float64)

    dydt[0] = -k1fwd * IL2Ra * IL2 + k1rev * IL2_IL2Ra - k6fwd * IL2Ra * IL2_gc + k6rev * IL2_IL2Ra_gc - k8fwd * IL2Ra * IL2_IL2Rb_gc + k8rev * IL2_IL2Ra_IL2Rb_gc - k12fwd * IL2Ra * IL2_IL2Rb + k12rev * IL2_IL2Ra_IL2Rb
    dydt[1] = -k2fwd * IL2Rb * IL2 + k2rev * IL2_IL2Rb - k7fwd * IL2Rb * IL2_gc + k7rev * IL2_IL2Rb_gc - k9fwd * IL2Rb * IL2_IL2Ra_gc + k9rev * IL2_IL2Ra_IL2Rb_gc - k11fwd * IL2Rb * IL2_IL2Ra + k11rev * IL2_IL2Ra_IL2Rb
    dydt[2] = -k3fwd * IL2 * gc + k3rev * IL2_gc - k5fwd * IL2_IL2Rb * gc + k5rev * IL2_IL2Rb_gc - k4fwd * IL2_IL2Ra * gc + k4rev * IL2_IL2Ra_gc - k10fwd * IL2_IL2Ra_IL2Rb * gc + k10rev * IL2_IL2Ra_IL2Rb_gc
    dydt[3] = -k11fwd * IL2_IL2Ra * IL2Rb + k11rev * IL2_IL2Ra_IL2Rb - k4fwd * IL2_IL2Ra * gc + k4rev * IL2_IL2Ra_gc + k1fwd * IL2 * IL2Ra - k1rev * IL2_IL2Ra
    dydt[4] = -k12fwd * IL2_IL2Rb * IL2Ra + k12rev * IL2_IL2Ra_IL2Rb - k5fwd * IL2_IL2Rb * gc + k5rev * IL2_IL2Rb_gc + k2fwd * IL2 * IL2Rb - k2rev * IL2_IL2Rb
    dydt[5] = -k6fwd *IL2_gc * IL2Ra + k6rev * IL2_IL2Ra_gc - k7fwd * IL2_gc * IL2Rb + k7rev * IL2_IL2Rb_gc + k3fwd * IL2 * gc - k3rev * IL2_gc
    dydt[6] = -k10fwd * IL2_IL2Ra_IL2Rb * gc + k10rev * IL2_IL2Ra_IL2Rb_gc + k11fwd * IL2_IL2Ra * IL2Rb - k11rev * IL2_IL2Ra_IL2Rb + k12fwd * IL2_IL2Rb * IL2Ra - k12rev * IL2_IL2Ra_IL2Rb
    dydt[7] = -k9fwd * IL2_IL2Ra_gc * IL2Rb + k9rev * IL2_IL2Ra_IL2Rb_gc + k4fwd * IL2_IL2Ra * gc - k4rev * IL2_IL2Ra_gc + k6fwd * IL2_gc * IL2Ra - k6rev * IL2_IL2Ra_gc
    dydt[8] = -k8fwd * IL2_IL2Rb_gc * IL2Ra + k8rev * IL2_IL2Ra_IL2Rb_gc + k5fwd * gc * IL2_IL2Rb - k5rev * IL2_IL2Rb_gc + k7fwd * IL2_gc * IL2Rb - k7rev * IL2_IL2Rb_gc
    dydt[9] = k8fwd * IL2_IL2Rb_gc * IL2Ra - k8rev * IL2_IL2Ra_IL2Rb_gc + k9fwd * IL2_IL2Ra_gc * IL2Rb - k9rev * IL2_IL2Ra_IL2Rb_gc + k10fwd * IL2_IL2Ra_IL2Rb * gc - k10rev * IL2_IL2Ra_IL2Rb_gc
    
    # added dydt[2] through dydt[9] based on the diagram pictured in type-I-ckine-model/model/graph.pdf on 9/19/17 by Adam; dydt[0] and dydt[1] were done by Aaron

    ###Ali's Addition
    



    return dydt

import numpy as np
from scipy.integrate import odeint

def dy_dt(y, t, IL2, IL9, k1fwd, k4fwd, k5rev, k6rev, k10rev, k11rev, k17rev, k18rev, k19rev):
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

    IL7Ra = y[10]
    IL7Ra_IL7 = y[11]
    gc_IL7 = y[12]
    IL7Ra_gc_IL7 = y[13]
    # k13 - k16

    IL9R = y[14]
    IL9R_IL9 = y[15]
    gc_IL9 = y[16]
    IL9R_gc_IL9 = y[17]
    # k17 - k20

    # The receptor-receptor forward rate is largely going to be determined by plasma membrane diffusion
    # so we're going to assume it's shared.
    k5fwd = k6fwd = k7fwd = k8fwd = k9fwd = k10fwd = k11fwd = k12fwd = k4fwd
    k19fwd = k20fwd = k4fwd

    # These are probably measured in the literature
    k1rev = k1fwd * 10 # doi:10.1016/j.jmb.2004.04.038, 10 nM
    k2fwd = k1fwd
    k2rev = k2fwd * 144 # doi:10.1016/j.jmb.2004.04.038, 144 nM
    k3fwd = k1fwd / 10.0 # Very weak, > 50 uM. Voss, et al (1993). PNAS. 90, 2428–2432.
    k3rev = 50000 * k3fwd

    k17fwd = k18fwd = k1fwd

    # To satisfy detailed balance these relationships should hold
    # _Based on initial assembly steps
    k4rev = k1fwd * k4fwd * k6rev * k3rev / k1rev / k6fwd / k3fwd
    k7rev = k3fwd * k7fwd * k2rev * k5rev / k2fwd / k5fwd / k3rev
    k12rev = k2fwd * k12fwd * k1rev * k11rev / k11fwd / k1fwd / k2rev
    # _Based on formation of full complex
    k9rev = k2rev * k10rev * k12rev / k2fwd / k10fwd / k12fwd / k3rev / k6rev * k3fwd * k6fwd * k9fwd
    k8rev = k2rev * k10rev * k12rev / k2fwd / k10fwd / k12fwd / k7rev / k3rev * k3fwd * k7fwd * k8fwd

    # _One detailed balance IL9 loop
    k20rev = k17rev * k19rev * k20fwd * k18fwd / k17fwd / k19fwd / k18rev

    dydt = np.zeros(y.shape, dtype = np.float64)

    dydt[0] = -k1fwd * IL2Ra * IL2 + k1rev * IL2_IL2Ra - k6fwd * IL2Ra * IL2_gc - k6rev * IL2_IL2Ra_gc - k8fwd * IL2Ra * IL2_IL2Rb_gc + k8rev * IL2_IL2Ra_IL2Rb_gc - k12fwd * IL2Ra * IL2_IL2Rb + k12rev * IL2_IL2Ra_IL2Rb
    dydt[1] = -k2fwd * IL2Rb * IL2 + k2rev * IL2_IL2Rb - k7fwd * IL2Rb * IL2_gc + k7rev * IL2_IL2Rb_gc - k9fwd * IL2Rb * IL2_IL2Ra_gc + k9rev * IL2_IL2Ra_IL2Rb_gc - k11fwd * IL2Rb * IL2_IL2Ra + k11rev * IL2_IL2Ra_IL2Rb
    dydt[2] = -k3fwd * IL2 * gc + k3rev * IL2_gc - k5fwd * IL2_IL2Rb * gc + k5rev * IL2_IL2Rb_gc - k4fwd * IL2_IL2Ra * gc + k4rev * IL2_IL2Ra_gc - k10fwd * IL2_IL2Ra_IL2Rb * gc + k10rev * IL2_IL2Ra_IL2Rb_gc
    dydt[3] = -k11fwd * IL2_IL2Ra * IL2Rb + k11rev * IL2_IL2Ra_IL2Rb - k4fwd * IL2_IL2Ra * gc + k4rev * IL2_IL2Ra_gc + k1fwd * IL2 * IL2Ra - k1rev * IL2_IL2Ra
    dydt[4] = -k12fwd * IL2_IL2Rb * IL2Ra + k12rev * IL2_IL2Ra_IL2Rb - k5fwd * IL2_IL2Rb * gc + k5rev * IL2_IL2Rb_gc + k2fwd * IL2 * IL2Rb - k2rev * IL2_IL2Rb
    dydt[5] = -k6fwd *IL2_gc * IL2Ra + k6rev * IL2_IL2Ra_gc - k7fwd * IL2_gc * IL2Rb + k7rev * IL2_IL2Rb_gc + k3fwd * IL2 * gc - k3rev * IL2_gc
    dydt[6] = -k10fwd * IL2_IL2Ra_IL2Rb * gc + k10rev * IL2_IL2Ra_IL2Rb_gc + k11fwd * IL2_IL2Ra * IL2Rb - k11rev * IL2_IL2Ra_IL2Rb + k12fwd * IL2_IL2Rb * IL2Ra - k12rev * IL2_IL2Ra_IL2Rb
    dydt[7] = -k9fwd * IL2_IL2Ra_gc * IL2Rb + k9rev * IL2_IL2Ra_IL2Rb_gc + k4fwd * IL2_IL2Ra * gc - k4rev * IL2_IL2Ra_gc + k6fwd * IL2_gc * IL2Ra - k6rev * IL2_IL2Ra_gc
    dydt[8] = -k8fwd * IL2_IL2Rb_gc * IL2Ra + k8rev * IL2_IL2Ra_IL2Rb_gc + k5fwd * gc * IL2_IL2Rb - k5rev * IL2_IL2Rb_gc + k7fwd * IL2_gc * IL2Rb - k7rev * IL2_IL2Rb_gc
    dydt[9] = k8fwd * IL2_IL2Rb_gc * IL2Ra - k8rev * IL2_IL2Ra_IL2Rb_gc + k9fwd * IL2_IL2Ra_gc * IL2Rb - k9rev * IL2_IL2Ra_IL2Rb_gc + k10fwd * IL2_IL2Ra_IL2Rb * gc - k10rev * IL2_IL2Ra_IL2Rb_gc
    
    dydt[10] = 0.0
    dydt[11] = 0.0
    dydt[12] = 0.0
    dydt[13] = 0.0

    dydt[2] = dydt[2] - k18fwd * IL9 * gc + k18rev * gc_IL9 - k19fwd * gc * IL9R_IL9 + k19rev * IL9R_gc_IL9
    dydt[14] = -k17fwd * IL9R * IL9 + k17rev * IL9R_IL9 - k20fwd * IL9R * gc_IL9 + k20rev * IL9R_gc_IL9
    dydt[15] = k17fwd * IL9R * IL9 - k17rev * IL9R_IL9 - k19fwd * gc * IL9R_IL9 + k19rev * IL9R_gc_IL9
    dydt[16] = -k20fwd * IL9R * gc_IL9 + k20rev * IL9R_gc_IL9 + k18fwd * IL9 * gc - k18rev * gc_IL9
    dydt[17] = k20fwd * IL9R * gc_IL9 - k20rev * IL9R_gc_IL9 + k19fwd * gc * IL9R_IL9 - k19rev * IL9R_gc_IL9

    # added dydt[2] through dydt[9] based on the diagram pictured in type-I-ckine-model/model/graph.pdf on 9/19/17 by Adam; dydt[0] and dydt[1] were done by Aaron

    return dydt

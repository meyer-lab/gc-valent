import numpy as np



def dy_dt(y, t, IL2, IL15, k4fwd, k5rev, k6rev, k10rev, k11rev):
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
    IL15Ra = y[10]
    IL15Rb = y[11]
    IL15_IL15Ra = y[12]
    IL15_IL15Rb = y[13]
    IL15_gc = y[14]
    IL15_IL15Ra_IL15Rb = y[15]
    IL15_IL15Ra_gc = y[16]
    IL15_IL15Rb_gc = y[17]
    IL15_IL15Ra_IL15Rb_gc = y[18]


    # The receptor-receptor forward rate is largely going to be determined by plasma membrane diffusion
    # so we're going to assume it's shared.
    k5fwd = k6fwd = k7fwd = k8fwd = k9fwd = k10fwd = k11fwd = k12fwd = k4fwd

    # These are probably measured in the literature
    k1fwd = 0
    k1rev = 0
    k2fwd = 0
    k2rev = 0
    k3fwd = 0
    k3rev = 0

    # To satisfy detailed balance these relationships should hold
    # _Based on initial assembly steps
    k4rev = k1fwd * k4fwd * k6rev * k3rev / k1rev / k6fwd / k3fwd
    k7rev = k3fwd * k7fwd * k2rev * k5rev / k2fwd / k5fwd / k3rev
    k12rev = k2fwd * k12fwd * k1rev * k11rev / k11fwd / k1fwd / k2rev
    # _Based on formation of full complex
    k9rev = k2rev * k10rev * k12rev / k2fwd / k10fwd / k12fwd / k3rev / k6rev * k3fwd * k6fwd * k9fwd
    k8rev = k2rev * k10rev * k12rev / k2fwd / k10fwd / k12fwd / k7rev / k3rev * k3fwd * k7fwd * k8fwd

    dydt = np.zeros(y.shape, dtype = np.float64)

    dydt[0] = -k1fwd * IL2Ra * IL2 + k1rev * IL2_IL2Ra - k6fwd * IL2Ra * IL2_gc - k6rev * IL2_IL2Ra_gc - k8fwd * IL2Ra * IL2_IL2Rb_gc + k8rev * IL2_IL2Ra_IL2Rb_gc - k12fwd * IL2Ra * IL2_IL2Rb + k12rev * IL2_IL2Ra_IL2Rb
    dydt[1] = -k2fwd * IL2Rb * IL2 + k2rev * IL2_IL2Rb - k7fwd * IL2Rb * IL2_gc + k7rev * IL2_IL2Rb_gc - k9fwd * IL2Rb * IL2_IL2Ra_gc + k9rev * IL2_IL2Ra_IL2Rb_gc - k11fwd * IL2Rb * IL2_IL2Ra + k11rev * IL2_IL2Ra_IL2Rb
    dydt[2] = 0
    dydt[3] = 0
    dydt[4] = 0
    dydt[5] = 0
    dydt[6] = 0
    dydt[7] = 0
    dydt[8] = 0
    dydt[9] = 0

    # TODO: Switch these to different rate constants
    dydt[10] = -k1fwd * IL15Ra * IL15 + k1rev * IL15_IL15Ra - k6fwd * IL15Ra * IL15_gc - k6rev * IL15_IL15Ra_gc - k8fwd * IL15Ra * IL15_IL15Rb_gc + k8rev * IL15_IL15Ra_IL15Rb_gc - k12fwd * IL15Ra * IL15_IL15Rb + k12rev * IL15_IL15Ra_IL15Rb
    dydt[11] = -k2fwd * IL15Rb * IL15 + k2rev * IL15_IL15Rb - k7fwd * IL15Rb * IL15_gc + k7rev * IL15_IL15Rb_gc - k9fwd * IL15Rb * IL15_IL15Ra_gc + k9rev * IL15_IL15Ra_IL15Rb_gc - k11fwd * IL15Rb * IL15_IL15Ra + k11rev * IL15_IL15Ra_IL15Rb
    dydt[12] = 0
    dydt[13] = 0
    dydt[14] = 0
    dydt[15] = 0
    dydt[16] = 0
    dydt[17] = 0
    dydt[18] = 0

    return dydt
import numpy as np



def dy_dt(y, t, IL2, k4fwd, k10rev):
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

    # The receptor-receptor forward rate is largely going to be determined by plasma membrane diffusion
    # so we're going to assume it's shared.
    k5fwd = k6fwd = k7fwd = k8fwd = k9fwd = k10fwd = k11fwd = k12fwd = k4fwd

    # To satisfy detailed balance these relationships should hold
    # TODO: Will be specified later
    # XXX: These assumptions may not be valid based on the literature
    k4rev = k5rev = k10rev # All involve gc coming off
    k6rev = k8rev = k12rev # All involve IL2Ra coming off
    k7rev = k9rev = k11rev # All involve IL2Rb coming off

    k11rev = k10rev * k2rev / k3rev # Same ratio of off rates as ligand
    k12rev = k10rev * k1rev / k3rev # Same ratio of off rates as ligand

    # These are measured in the literature
    k1fwd = 0
    k1rev = 0
    k2fwd = 0
    k2rev = 0
    k3fwd = 0
    k3rev = 0


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

    return dydt
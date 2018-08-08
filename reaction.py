import sys
import math

import const as cs

###############################################################################
# Define reaction path 
###############################################################################
class Reaction(object):

###############################################################################
# name - name of reaction, i.e dissociation, transition, rotation, etc
# in_state/ts_state/fi_state - configurations along reaction coordinate
# T_eV - temperature in eV
###############################################################################
    def __init__(self, name, in_state, ts_state, fi_state, T_K):
        self.name = name
        self.in_state = in_state
        self.ts_state = ts_state
        self.fi_state = fi_state
        self.T_eV = cs.k_Bol * T_K 
        return None

###############################################################################
# Calculate DFT activation energy of the reaction 
###############################################################################
    def calc_E_DFT_activ(self):
        E_DFT_activ = self.ts_state.E_DFT - self.in_state.E_DFT
        return E_DFT_activ

###############################################################################
# Calculate Gibbs activation energy at T = 0 K
###############################################################################
    def calc_Gibbs_0_activ(self):
        Gibbs_0_activ = self.ts_state.Gibbs_0 - self.in_state.Gibbs_0
        return Gibbs_0_activ

###############################################################################
# Calculate Gibbs activation energy at finite temperature
###############################################################################
    def calc_Gibbs_T_activ(self, T_K):
        Gibbs_T_activ = (self.ts_state.calc_Gibbs_T(T_K) 
                         - self.in_state.calc_Gibbs_T(T_K))
        return Gibbs_T_activ

###############################################################################
# Calculate DFT reaction energy
###############################################################################
    def calc_E_DFT_react(self):
        E_DFT_react = self.fi_state.E_DFT - self.in_state.E_DFT
        return E_DFT_react

###############################################################################
# Calculate Gibbs reaction free energy at zero temperature
###############################################################################
    def calc_Gibbs_0_react(self):
        Gibbs_0_react = self.fi_state.Gibbs_0 - self.in_state.Gibbs_0
        return Gibbs_0_react

###############################################################################
# Calculate Gibbs reaction free energy at finite temperature
###############################################################################
    def calc_Gibbs_T_react(self, T_K):
        Gibbs_T_react = (self.fi_state.calc_Gibbs_T(T_K) 
                         - self.in_state.calc_Gibbs_T(T_K))
        return Gibbs_T_react

###############################################################################
# Calculate tunneling coefficient
###############################################################################
    def calc_tunnel_coef(self, T_K):
        try:
            freq_im_max = max(self.ts_state.freq_im)
        except ValueError:
            sys.exit("The transition state has no imaginary frequencies")
        beta = 2.0 * math.pi * self.calc_E_DFT_activ() / freq_im_max
        Upp = freq_im_max / self.T_eV
        k1 = 0.5 * Upp / math.sin(Upp / 2.0)
        k2 = (Upp * math.exp(-self.calc_E_DFT_activ() 
                             / self.T_eV ) * math.exp(-beta)) 
        k3 = (1.0 / (2.0 * math.pi - Upp) 
              - math.exp(-beta) / (4.0 * math.pi - Upp))
        kappa =  k1 - k2 * k3
        return kappa

###############################################################################
# Calculate rate of the reaction
###############################################################################
    def calc_rate_const(self, T_K):
        rate_const = (self.calc_tunnel_coef(T_K) * self.T_eV 
                      * math.exp(-self.calc_Gibbs_T_activ(T_K) / self.T_eV) 
                      / cs.h_bar)
        return rate_const

###############################################################################
# Print details of reaction and calculate missing values
###############################################################################
    def __str__(self):
        T_K = self.T_eV / cs.k_Bol
        str_ = ("{:-^66s}\n"
                "Termperature:                                  {:>16.0f}   K\n"
                "Finite-temperature Gibbs free energies:\n"
                "Initial                                        {:>16.3f}  eV\n"
                "Transition                                     {:>16.3f}  eV\n"
                "Final                                          {:>16.3f}  eV\n"
                "DFT activation energy:                         {:>16.3f}  eV\n"
                "Zero-temperature activation Gibbs energy:      {:>16.3f}  eV\n"
                "Finite-temperature activation Gibbs energy:    {:>16.3f}  eV\n"
                "DFT reaction energy:                           {:>16.3f}  eV\n"
                "Zero-temperature reaction Gibbs energy:        {:>16.3f}  eV\n"
                "Finite-temperature reaction Gibbs energy:      {:>16.3f}  eV\n"
                "Tunneling factor:                              {:>16.3f}\n"
                "Rate constant:                                 {:>16.0e} 1/s\n"
                "------------------------------------------------------------"
                "------\n"
                .format(self.name, 
                        T_K, 
                        self.in_state.calc_Gibbs_T(T_K), 
                        self.ts_state.calc_Gibbs_T(T_K), 
                        self.fi_state.calc_Gibbs_T(T_K),
                        self.calc_E_DFT_activ(),
                        self.calc_Gibbs_0_activ(),
                        self.calc_Gibbs_T_activ(T_K),
                        self.calc_E_DFT_react(),
                        self.calc_Gibbs_0_react(),
                        self.calc_Gibbs_T_react(T_K),
                        self.calc_tunnel_coef(T_K),
                        self.calc_rate_const(T_K)) )
        return str_



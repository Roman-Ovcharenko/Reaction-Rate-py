import math

import const as cs

###############################################################################
# Define configuration according to the transition state theory
# see H. Eyring, J. Chem. Phys. 3, 107 (1935).
###############################################################################
class Config(object):

###############################################################################
# type = initial/final/transition state
# file_name - where data was taken from 
# E_DFT - the total DFT electronic energy of the system at T = 0 K
# frec_re - real harmonic vibrational frequencies
# frec_im - imaginary harmonic vibrational frequencies along reaction
# coordinate
###############################################################################
    def __init__(self, typ):
        self.type = typ
        self.file_name = self._get_file_name()
        self.E_DFT, self.freq_re, self.freq_im = self._get_param()
        self.Gibbs_0 = self._calc_Gibbs_0()
        return None

###############################################################################
# read name of the file for current configuration
###############################################################################
    def _get_file_name(self):
        return input("Type file name for the {} configuration: "
                     .format(self.type))

###############################################################################
# Read the total DFT energy and real/imaginary vibrational frequencies 
# from the input file
###############################################################################
    def _get_param(self):
        freq_re = []
        freq_im = []
        with open(self.file_name,'r') as f:
            for line in f:
                if "free  energy   TOTEN  =" in line:
                    lst = line.split()
                    E_DFT = float(lst[4])
                if "f  =" in line:
                    lst = line.split()
                    freq_re.append(float(lst[9]) * cs.meV_to_eV)
                if "f/i= " in line:
                    lst = line.split()
                    freq_im.append(float(lst[8]) * cs.meV_to_eV)
        return E_DFT, freq_re, freq_im

###############################################################################
# Print an instance of the configuration object 
###############################################################################
    def __str__(self):
        str_ = ("{:-^55s}\n"
                "Related file:                      {:>20s}\n"
                "Total DFT energy:                   {:>16.3f} eV\n"
                "Gibbs zero temperature free energy: {:>16.3f} eV\n"
                "-------------------------------------------------------\n"
                .format(self.type, self.file_name, self.E_DFT, self.Gibbs_0))
        return str_

###############################################################################
# Zero temperature correction to DFT energy
###############################################################################
    def _calc_Gibbs_0(self):
        del_Gibbs_0 = math.fsum(self.freq_re) / 2.0
        Gibbs_0 = self.E_DFT + del_Gibbs_0
        return Gibbs_0

###############################################################################
# Finite temperature correction to Gibbs free energy at T = 0 K 
###############################################################################
    def calc_Gibbs_T(self, T_K):
        T_eV = cs.k_Bol * T_K
        del_Gibbs_T = 0.0
        for nu in self.freq_re:
            del_Gibbs_T = (del_Gibbs_T 
                           + T_eV * math.log(1.0 - math.exp(-nu/T_eV)))
        Gibbs_T = self.Gibbs_0 + del_Gibbs_T 
        return Gibbs_T


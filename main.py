#!/usr/local/bin/python3.7

import const as cs
from configuration import Config
from reaction import Reaction

# Read temperature
print("Calculation of the reaction rate")
T_K = float(input("Type temperature (in K): "))
if T_K < 0.0:
    raise Exception("Negative temperature")

# Initiate thermodynamics objects
in_state = Config("initial state")
ts_state = Config("transition state")
fi_state = Config("final state")
print()

print(in_state)
print(ts_state)
print(fi_state)

diss = Reaction("dissociation", in_state, ts_state, fi_state, T_K)
print(diss)
       



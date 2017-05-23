import numpy as np
from Model import InitialConditionsUtils as ICU
from Initial import Computation as CPT

class InitialConditions(object):

    def __init__(self, percentage_loaded,  soil_depth, f, BC, SD, HS, SO, Init=0):

      self.recharge_rate = CPT.compute_recharge_rate(SO.t[0], SO.t, SO.recharge_type, SO.period, SO.recharge_chronicle)
      self.Smax = ICU.compute_Smax(f, SD, soil_depth)
      self.sin = ICU.compute_sin(percentage_loaded, self.Smax, BC, Init)
      self.qin = ICU.compute_qin(self.sin, f, k , BC, SD, HS, Init)
      self.q_sin = ICU.compute_q_sin(self.sin, self.qin, SO.t, f, SD, HS, self.recharge_rate, Init)

    def get_sin():
        return self.sin

    def get_qin():
        return self.qin

    def get_q_sin():
        return self.q_sin

    def get_Smax():
        return self.Smax

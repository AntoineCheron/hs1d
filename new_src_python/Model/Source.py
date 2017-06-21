from Model import TimePropertiesUtils as TP
from Model import TimeUnit as TU
import numpy as np

class Source(object):
  """
      Class managing temporal aspect and source terms (recharge)

      @param
        period : period used to compute the recharge chronicle (if not 'databased')
        recharge_type : type of used recharge ('periodical, 'square', 'steady', 'random'
                                               'databased')
        recharge_rate : either a float or a vector of recharge (m/s)
        tmin : minimal time value of the serie
        tmax : maximal time value of the serie
        Nt : number of time values
        unit : time unit ('days', 'hour', 'min', 'sec', 'year')
        time_custom : -1 or a vector containing time values if 'databased' recharge

      @attributes
        recharge : recharge value at the asked time
        period : period used to compute the recharge chronicle (if not 'databased')
        tmax : maximal time value of the chronicle
        recharge_type : type of used recharge
        t_chronicle : a vector containing time values of the chronicle
        recharge_chronicle : all recharge values
        TP : TimeProperties class

  """

  def __init__(self, period, recharge_type, recharge_rate, tmin, tmax, Nt, unit, time_custom):
      self.recharge_type = recharge_type
      self.t, self.tmax, self.tmin = TP.time_properties(tmin, tmax, Nt, unit, time_custom)
      self.period = self.source_terms(period, recharge_type)
      self.recharge_chronicle = self.set_recharge_chronicle(recharge_rate, self.period, self.t, recharge_type)
    
  def source_terms(self, period, recharge_type):
    
      if recharge_type == 'periodical' or recharge_type == 'square':
          if period is None:
              return TU.time_to_seconds(5)

      elif recharge_type == 'steady' or recharge_type is None:
          return 'inf'
    
      elif recharge_type == 'random':
          return 0
    
      elif recharge_type == 'databased':
          return -1
    
      return self.period
    
  def set_recharge_chronicle(self, recharge_rate, period, t, recharge_type):
    
      #        print('time=',self.TP.t)
      if period == 'inf' or period is None:
          recharge_rate = (recharge_rate*10**-3)/86400
          return recharge_rate*np.ones(t)
    
      elif period == 0:
          return self.set_random_recharge()
    
      elif period == -1:
          return recharge_rate
    
      else:
          recharge_rate = (recharge_rate*10**-3)/86400
          if recharge_type == 'periodical':
              return recharge_rate*(1+np.cos(2*np.pi*(t/period)))
          elif recharge_type == 'square':
              Int_ = np.floor(t/period)
              Odd_rest = Int_%2
              recharge_chronicle = recharge_rate*Odd_rest
              rem_stiff = t%(2*period)
              bool1 = rem_stiff < 60
              recharge_chronicle[bool1] = recharge_rate*(1-(rem_stiff[bool1])/60)
              rem_stiff = (t + period)%(2*period)
              bool1 = rem_stiff < 60
              recharge_chronicle[bool1] = recharge_rate*((rem_stiff[bool1])/60)
              return recharge_chronicle
    
    ###########################################################################
    #                               GETTERS                                   #
    ###########################################################################

  def get_t(self):
        return self.t
    
  def get_recharge(self):
        return self.recharge_chronicle
    
  def get_recharge_type(self):
        return self.recharge_type
    
  def get_period(self):
        return self.period

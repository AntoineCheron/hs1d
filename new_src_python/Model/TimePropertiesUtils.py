import numpy as np
import TimeUnit as TU

"""
    Properties management for time (unit and values) used to describe the
    hillslope and its paramters variations
    #######################################################################
    @param
        tmin : minimal time value (usually 0)
        tmax : maximal time value
        Nt : number of time steps from tmin to tmax
        unit : string corresponding to the time unit of the time serie
            choices : year, days, hour, min, sec
        time_custom : -1 or a vector of time values if databased recharge

    @attributes
      t : time vector
      Nt : number of time values (=len(t))
      tmax : maximal time value
      tmin : minimal time value
      unit : time unit
      TU : TimeUnit class
    #######################################################################
"""

def time_properties(tmin=0, tmax=35, Nt=35*24*10, unit='days', time_custom=-1):
    """
        Creates a time vector in sec, based on tmin, tmax and Nt (after
        conversion from unit to seconds)
    """
    TU = TU.TimeUnit(tmax, unit)

    # Convert tmax and tmin to seconds
    tmax = TU.time_to_seconds(tmax)
    tmin = TU.time_to_seconds(tmin)

    if isinstance(time_custom, int):
        t = np.arange(tmin, tmax, (tmax-tmin)/(Nt-1))
        t = np.reshape(t,(len(t),1))
    else :
        t = TU.time_to_seconds(time_custom)
    return t

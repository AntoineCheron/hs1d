"""
    Module managing the conversion between time units dor the calculation
    ######################################################################
    @param
      tmax : the maximal time of the timeserie
      unit : string corresponding to the time unit
        choices : year, days, hour, min, sec

    @attributes
      tmax : the maximal time of the timeserie
      unit : string corresponding to the time unit
    ######################################################################
"""

def time_to_seconds(tmax, unit):
    """
        Conversion between each time unit to seconds
    """
    sec_in_min = 60
    sec_in_hour = 3600
    sec_in_day = 86400
    sec_in_year = 86400*365

    if unit == '':
        unit = ('enter time units - seconds, minutes, days, year')
        if unit == 'days':
            return tmax*sec_in_day
        elif unit == 'hour':
            return tmax*sec_in_hour
        elif unit == 'min':
            return tmax*sec_in_min
        elif unit == 'sec':
            return tmax
        elif unit == 'year':
            return tmax*sec_in_year
        else:
            print('UNIT UNKNOWN #s\n', unit)
            return tmax

    elif unit == 'year':
        return tmax*sec_in_year
    elif unit == 'days':
        return tmax * sec_in_day
    elif unit == 'hour':
        return tmax * sec_in_hour
    elif unit == 'min':
        return tmax * sec_in_min
    elif unit == 'sec':
        return tmax
    else:
        print('UNIT UNKNOWN #s\n', unit)
        return tmax

def time_to_days(tmax, unit):
    sec_in_day = 86400
    min_in_day = 1440
    hour_in_day = 24
    day_in_year = 365


    if unit == '':
        unit = ('enter time units - seconds, minutes, days, year')
        if unit == 'hour':
            return tmax/hour_in_day
        elif unit == 'min':
            return tmax/min_in_day
        elif unit == 'sec':
            return tmax/sec_in_day
        elif unit == 'year':
            return tmax*day_in_year
        else:
            print('UNIT UNKNOWN #s\n', unit)
            return tmax

    elif unit == 'sec':
        return tmax/sec_in_day
    elif unit == 'min':
        return tmax / min_in_day
    elif unit == 'hour':
        return tmax / hour_in_day
    elif unit == 'year':
        return tmax*day_in_year
    else:
        print('UNIT UNKNOWN #s\n', unit)
        return tmax

def time_to_years(tmax, unit):
    days_in_year = 365
    hours_in_year = 24*days_in_year
    min_in_year = 60*hours_in_year
    sec_in_year = 60*min_in_year


    if unit == '':
        unit = ('enter time units - seconds, minutes, days, year')
        if unit == 'days':
            return tmax*days_in_year
        elif unit == 'hour':
            return tmax*hours_in_year
        elif unit == 'min':
            return tmax*min_in_year
        elif unit == 'sec':
            return tmax
        else:
            print('UNIT UNKNOWN #s\n', unit)
            return tmax

    elif unit == 'sec':
        return tmax / sec_in_year
    elif unit == 'min':
        return tmax / min_in_year
    elif unit == 'hour':
        return tmax / hours_in_year
    elif unit == 'days':
        return tmax / days_in_year
    else:
        print('UNIT UNKNOWN #s\n', unit)
        return tmax

def time_to_hours(tmax, unit):
    hours_in_day = 24
    hours_in_year = 365*hours_in_day
    min_in_hour = 60
    sec_in_hour = 60*min_in_hour

    if unit == '':
        unit = ('enter time units - seconds, minutes, days, year')
        if unit == 'year':
            return tmax*hours_in_year
        elif unit == 'days':
            return tmax*hours_in_day
        elif unit == 'hour':
            return tmax
        elif unit == 'min':
            return tmax/min_in_hour
        elif unit == 'sec':
            return tmax/sec_in_hour
        else:
            print('UNIT UNKNOWN #s\n', unit)
            return tmax

    elif unit == 'sec':
        return tmax / sec_in_hour
    elif unit == 'min':
        return tmax / min_in_hour
    elif unit == 'hour':
        return tmax
    elif unit == 'days':
        return tmax * hours_in_day
    elif unit == 'year':
        return tmax * hours_in_year
    else:
        print('UNIT UNKNOWN #s\n', unit)
        return tmax

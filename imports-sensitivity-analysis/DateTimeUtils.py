"""
Created on Wed Sep 16 17:33:58 2020

Functions for converting between calandar dates and decimal (float) times

@author: david
"""
from dateutil import parser
from datetime import datetime, timedelta
from time import mktime

def date2FloatYear(dates):
    
    def sinceEpoch(date): # returns seconds since epoch
        return mktime(date.timetuple())
    s = sinceEpoch

    "Convert to list if not in list"
    if not isinstance(dates, list):    
        dates = [dates]
    
    "Convert to datetime objects if not datetimes"
    if not isinstance(dates[0], datetime):
        dates = [parser.parse(dt) for dt in dates]
    
    for idx, date in enumerate(dates):
        year = date.year
        startOfThisYear = datetime(year=year, month=1, day=1)
        startOfNextYear = datetime(year=year+1, month=1, day=1)
        yearElapsed = s(date) - s(startOfThisYear)
        yearDuration = s(startOfNextYear) - s(startOfThisYear)
        fraction = yearElapsed/yearDuration
        dates[idx] = date.year + fraction
        
    if len(dates) == 1:
        dates = dates[0]

    return dates

def floatYear2Date(times):
    
    if not isinstance(times, list):    
        times = [times]
    
    for idx, time in enumerate(times):
        year = int(time)
        rem = time - year
        base = datetime(year, 1, 1)
        time = base + timedelta(seconds=(base.replace(year=base.year + 1) - base).total_seconds() * rem)    
        times[idx] = time.date()
        
    if len(times) == 1:
        times = times[0]
        
    return times

if  __name__ == '__main__':
    
    "Test conversion to float and back to datetime object"
    date = parser.parse("2020-01-01")
    print('Date: ' + str(date))
    
    float_year = date2FloatYear(date)
    print('Float year: ' + str(float_year))
    
    date = floatYear2Date(float_year)
    print('Converted date: ' + str(date))
    
    "Test with lists"
    date1 = parser.parse("2020-01-01")
    date2 = parser.parse("2020-03-19")
    dates = [date1,date2]
    
    float_years = date2FloatYear(dates)
    print('Float years: ' + str(float_years[0]) + ' ' + str(float_years[1]))
    
    dates = floatYear2Date(float_years)
    print('Converted dates: ' + str(dates[0]) + ' ' + str(dates[1]))
    

    
    
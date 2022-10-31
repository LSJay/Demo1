import datetime
import sys
import os
import subprocess
from os.path import join as path
file = sys.argv[1]
'''
def my_func(year = 2015,month = 1,day = 23):
    atime = datetime.date(year, month, day)
    b = datetime.time(12,00,00)
    print(atime,b)



c = list(sys.argv[1:])
#L = ['year','month','day']
year = int(c[0])
print(year)
month = int(c[1])
print(month)
day = int(c[2])
print(day)


my_func(year,month,day)
'''


commod =  ' awk -F "=" "{printf $5}" file'
subprocess.call(commod,shell=True)

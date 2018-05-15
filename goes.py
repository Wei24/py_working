import os 
from sunpy.time import TimeRange 
from sunpy import lightcurve as lc 
import pickle 
goesplottim = TimeRange("2017/09/06 18:59:20.000", "2017/09/06 19:59:28.000") 
goes = lc.GOESLightCurve.create(goesplottim) 
fi2 = open("./goes.dat", "wb") 
pickle.dump(goes, fi2) 
fi2.close()
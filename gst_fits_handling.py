#trying to co-align gst data with the NORMAL ones :-\
def de_saturat()
    [rows,cols]=fimage.shape
    for i in range(rows-1):
        for j in range(cols-1):
            if fimage[j,i] < 0:
                
         
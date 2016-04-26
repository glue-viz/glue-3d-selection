import time  
import Image
import flood_fill
import numpy as np
# a useful website for using the .pyx and test the time is http://blog.perrygeo.net/2008/04/19/a-quick-cython-introduction/
# the timeit module could be use for tracing the running time of the code http://www.diveintopython.net/performance_tuning/timeit.html

jpgfile = Image.open('zoo_test.jpg').convert('L') # convert RGB img to greyscale
data = np.array(jpgfile)

row = data.shape[0]
col = row = data.shape[1]

newdata = []
start = time.time() # add time tracker
newdata =flood_fill.cyfill(data, (row/2, col/2), 255, 30.0) # (data, start_coords,fill_value, threshold)
# TODO: threshold here ranges from greyscale 0~255 but should be converted to be related with mouse dragging distance

end = time.time()
newdata = np.array(newdata)
newim = Image.fromarray(newdata)
newim.save('new_zoo_test', 'JPEG')
print('New image produced within '+str(end-start)+' sec')

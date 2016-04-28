import time  
import Image
import flood_fill_3d
import numpy as np
from astropy.io import fits
# a useful website for using the .pyx and test the time is http://blog.perrygeo.net/2008/04/19/a-quick-cython-introduction/
# the timeit module could be use for tracing the running time of the code http://www.diveintopython.net/performance_tuning/timeit.html

hdulist = fits.open('l1448_13co.fits')
fitsdata = hdulist[0].data
# fitsdata.dtype = '>f4', which can't be transferred into Cython as input :( So we use a test data

x = int(fitsdata.shape[0])
y = int(fitsdata.shape[1])
z = int(fitsdata.shape[2])

# create a test data here
# data = np.rint(np.random.rand(x, y, z)*40)
# data[x//2-10:x//2+10, y//2-10:y//2+10, z//2-10:z//2+10] = 0
data = np.asarray(fitsdata, np.float64)

oldhdu = fits.PrimaryHDU(data)
oldhdulist = fits.HDUList([oldhdu])
oldhdulist.writeto('old.fits')

newdata = []
start = time.time() # add time tracker
newdata =flood_fill_3d.cyfill(data, (int(x//2), int(y//2), int(z//2)), 5, 1) # (data, start_coords,fill_value, threshold)
# TODO: threshold here ranges from greyscale 0~255 but should be converted to be related with mouse dragging distance

end = time.time()
newdata = np.array(newdata)

# create a new image file
newhdu = fits.PrimaryHDU(newdata)
newhdulist = fits.HDUList([newhdu])
newhdulist.writeto('new.fits')

print('New image produced within '+str(end-start)+' sec')

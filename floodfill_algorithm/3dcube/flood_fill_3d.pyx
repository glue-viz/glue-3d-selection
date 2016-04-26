import numpy as np
cimport numpy as np
# def cyfill(float[:, :, :] data,   #np.dtype('>f4')
def cyfill(np.ndarray[np.float64_t, ndim=3] data, tuple start_coords, Py_ssize_t fill_value, float threshold):
    """
    Flood fill algorithm
    
    Parameters
    ----------
    data : (M, N) ndarray of uint8 type
        Image with flood to be filled. Modified inplace.
    start_coords : tuple
        Length-2 tuple of ints defining (row, col) start coordinates.
    fill_value : int
        Value the flooded area will take after the fill.
        
    Returns
    -------
    None, ``data`` is modified inplace.
    """
    cdef:
        Py_ssize_t x, y, z, xsize, ysize, zsize, orig_value
        float thres 
        set stack
    
    xsize = data.shape[0]
    ysize = data.shape[1]
    zsize = data.shape[2]
    
    print data.shape[0]
    print data[20, 20, 20]
    orig_value = int(data[start_coords[0], start_coords[1], start_coords[2]])
    # orig_value = data[20, 20, 20]
	# this threshold should combine the relation between dragging distance and greyscale(0~255)
    thres = threshold
    if fill_value == orig_value:
        raise ValueError("Filling region with same value "
                         "already present is unsupported. "
                         "Did you already fill this region?")
    
    stack = set(((start_coords[0], start_coords[1], start_coords[2]),))

    while stack:
        x, y, z = stack.pop()
		# set the threshold as 30 here but this could be get through an input
        if data[x, y, z]>orig_value-thres and data[x, y, z]<orig_value+thres:
            data[x, y, z] = fill_value
            if x > 0:
                stack.add((x - 1, y, z))
            if x < (xsize - 1):
                stack.add((x + 1, y, z))
            if y > 0:
                stack.add((x, y - 1, z))
            if y < (ysize - 1):
                stack.add((x, y + 1, z))
            if z > 0:
                stack.add((x, y, z - 1))
            if z < (zsize - 1):
                stack.add((x, y, z + 1))
    return data
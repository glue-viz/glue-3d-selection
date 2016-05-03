from astropy.io import fits
from floodfill_scipy import floodfill_scipy
import matplotlib.pyplot as plt
from flood_fill_3d import cyfill

data = fits.open('slice.fits')[0].data.astype(float)
start = (0, 57, 11)

for algorithm in ['penny', 'tom']:

    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(3, 3, 1)
    ax.imshow(data[0], cmap=plt.cm.viridis, vmin=0, vmax=4, origin='lower',interpolation='nearest')
    ax.set_autoscale_on(False)
    ax.plot(start[2], start[1], 'ro', markersize=5)
    ax.xaxis.set_ticklabels('')
    ax.yaxis.set_ticklabels('')

    for i, threshold in enumerate([1.001, 1.01, 1.05, 1.1, 1.2, 1.5, 1.76, 2]):

        if algorithm == 'tom':
            mask = floodfill_scipy(data, start, threshold)
        else:
            mask = data.copy()
            cyfill(mask, start, 1, threshold)
            mask = mask == 1

        ax = fig.add_subplot(3, 3, 2 + i)
        ax.text(0.05, 0.95, str(threshold), size=10, color='white', transform=ax.transAxes, ha='left', va='top')
        ax.imshow(mask[0], cmap=plt.cm.gist_heat, vmin=0, vmax=2, origin='lower',interpolation='nearest')

        ax.xaxis.set_ticklabels('')
        ax.yaxis.set_ticklabels('')

    fig.savefig('test_case_{0}.png'.format(algorithm), dpi=150)

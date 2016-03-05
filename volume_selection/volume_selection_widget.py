"""
    3D Volume Selection tool based on Vispy library
    Date: Feb 15, 2016

    Usage: running the code with 'python selection_vispy_github.py'.

    Keypress:  1 for free lasso selection
               2 for rectangle selection
               3 for ellipse selection
               4 for individual point picking
"""

import sys
import numpy as np

from PyQt4 import QtGui, QtCore
from vispy import app, scene
from volume_selection_common import SelectionCommon
import pyfits

from multivol import MultiVolume
from multivol import get_translucent_cmap

app.use_app('pyqt4')


class DemoScene(QtGui.QWidget):
    def __init__(self, keys='interactive'):
        super(DemoScene, self).__init__()
        # Layout and canvas creation
        box = QtGui.QVBoxLayout(self)
        self.resize(500,500)
        self.setLayout(box)
        self.canvas = scene.SceneCanvas(keys=keys)
        box.addWidget(self.canvas.native)

        # Camera
        self.view = self.canvas.central_widget.add_view()
        self.view.camera = scene.cameras.TurntableCamera(elevation = 25, azimuth=20, distance = 2.0, center=(0,0,0))

        # Data
        fitsdata = pyfits.open('l1448_13co.fits')
        vol1 = np.nan_to_num(fitsdata[0].data)
        self.vol_data = vol1

        """
        The transpose here and after is for solving the coordinate mismatch between volume visual input data and its
        rendering result. The rendered volume shown on 2D screen, or 'what we see', is through displaying transform
        (self.tr here) of 'transposed input data', thus we use 'transpose' to enable our selection focusing on 'what
        we see' on the screen rather than the real input data of volume.

        """
        new_pos = np.transpose(vol1)

        # TODO: replace the min&max threshold with real settings in Glue UI
        min_threshold = np.min(self.vol_data)
        max_threshold = np.max(self.vol_data)
        self.pos_data = np.argwhere(new_pos >= min_threshold)  # get voxel positions

        grays = get_translucent_cmap(1, 1, 1)

        self.volume_pool = [(vol1, (1, 6), grays)]
        self.volume = MultiVolume(self.volume_pool)
        self.trans = [-vol1.shape[2]/2., -vol1.shape[1]/2., -vol1.shape[0]/2.]
        self.volume.transform = scene.STTransform(translate=self.trans)
        self.view.add(self.volume)

        # Add a 3D axis to keep us oriented
        axis = scene.visuals.XYZAxis(parent=self.view.scene)

        self.selection = SelectionCommon(canvas=self.canvas, view=self.view,
                                          vol_data=self.vol_data, volume=self.volume,
                                          volume_pool=self.volume_pool, pos_data=self.pos_data)

if __name__ == '__main__':
    appQt = QtGui.QApplication(sys.argv)
    view = DemoScene(keys='interactive')
    view.show()
    appQt.exec_()

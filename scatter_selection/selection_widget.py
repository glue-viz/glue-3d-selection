"""
    3D Points Picking & Selection tool based on Vispy library
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
from selection_common import SelectionCommon

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

        # Test data
        self.data = np.random.uniform(-1, 1, size=(200, 3))
        # TODO: set facecolor and ptsize from Glue UI control
        self.facecolor = np.ones((200,4), dtype=np.float)
        self.ptsize = 10

        # Draw test points data
        self.scatter = scene.visuals.Markers()
        self.scatter.set_data(self.data, face_color=self.facecolor, edge_color=None, size=self.ptsize)
        self.view.add(self.scatter)

        # Add a 3D axis to keep us oriented
        axis = scene.visuals.XYZAxis(parent=self.view.scene)

        self.selection = SelectionCommon(canvas=self.canvas, view=self.view,
                                          data=self.data, scatter=self.scatter,
                                          facecolor=self.facecolor)

if __name__ == '__main__':
    appQt = QtGui.QApplication(sys.argv)
    view = DemoScene(keys='interactive')
    view.show()
    appQt.exec_()

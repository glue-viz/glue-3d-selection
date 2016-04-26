__author__ = 'penny'

"""
    3D Ray Cursor based on Vispy library
    Date: May 6, 2016
    Usage:
    Keypress '1': 'ray_cursor'
    Press again to switch between selection & view mode :)

"""

import sys
import numpy as np

from PyQt4 import QtGui, QtCore
from vispy import app, scene, visuals
from astropy.io import fits

app.use_app('pyqt4')

class DemoScene(QtGui.QWidget):
    def __init__(self, keys='interactive'):
        super(DemoScene, self).__init__()

        # Layout and canvas creation
        box = QtGui.QVBoxLayout(self)
        self.resize(800, 600)
        self.setLayout(box)
        self.canvas = scene.SceneCanvas(keys=keys)
        box.addWidget(self.canvas.native)

        # Connect events
        self.canvas.events.mouse_press.connect(self.on_mouse_press)
        self.canvas.events.key_press.connect(self.on_key_press)

        # Setup some defaults
        self.mesh = None
        self.selected = []
        self.white = (1.0, 1.0, 1.0, 1.0)
        self.black = (0.0, 0.0, 0.0, 0.0)

        # Camera
        self.view = self.canvas.central_widget.add_view()
        self.view.camera = scene.cameras.TurntableCamera(elevation = 90, azimuth=0, fov=60, center=(0,0,0), distance=10.)

        # Data
        fitsdata = fits.open('test_data/l1448_13co.fits')
        self.vol_data = np.nan_to_num(fitsdata[0].data)

        new_pos = np.transpose(self.vol_data)

        # TODO: replace the min&max threshold with real settings in Glue UI
        self.pos_data = np.indices(self.vol_data.shape).reshape(3,-1).transpose()

        self.volume = scene.visuals.Volume(self.vol_data, parent=self.view.scene)
        self.trans = [-self.vol_data.shape[2]/2., -self.vol_data.shape[1]/2., -self.vol_data.shape[0]/2.]
        self.volume.transform = scene.STTransform(translate=self.trans)

        self.tr = self.volume.get_transform(map_from='visual', map_to='canvas')
        self.tr_back = self.volume.get_transform(map_from='canvas', map_to='visual')

        # Add a text instruction
        self.text = scene.visuals.Text('', color='red', pos=(self.canvas.size[0]/4.0,  20), parent=self.canvas.scene)

        # Add a cursor line
        self.ray_line = scene.visuals.Line(color='green', width=10, parent=self.view.scene)
        # Add text instruction for start point and end point
        self.rayline_text = scene.visuals.Text('', color='green', font_size=10, pos=(self.canvas.size[0]/2.0, 500),
                                               parent=self.canvas.scene)

        # Add a 3D axis to keep us oriented
        axis = scene.visuals.XYZAxis(parent=self.view.scene)

        # Selection
        self.selection_flag = False
        self.selection_pool = {'1': 'ray_cursor'}
        self.selection_id = '1'  # default as 1
        self.selection_origin = (0, 0)

    def transform(self, data):  # from visual to canvas
        data = self.tr.map(data)
        data /= data[:, 3:]  # normalize with homogeneous coordinates
        return data[:, :2]

    def transform_back(self, data):  # from canvas to visual
        data = self.tr_back.map(data)
        data = data[:3]/data[3]  # normalize with homogeneous coordinates
        return data

    def event_connect(self, flag):
        if flag:
            self.view.camera._viewbox.events.mouse_move.disconnect(
                    self.view.camera.viewbox_mouse_event)
            self.view.camera._viewbox.events.mouse_press.disconnect(
                        self.view.camera.viewbox_mouse_event)
            self.view.camera._viewbox.events.mouse_release.disconnect(
                        self.view.camera.viewbox_mouse_event)
            self.view.camera._viewbox.events.mouse_wheel.disconnect(
                        self.view.camera.viewbox_mouse_event)
        else:
            self.view.camera._viewbox.events.mouse_move.connect(
                    self.view.camera.viewbox_mouse_event)
            self.view.camera._viewbox.events.mouse_press.connect(
                        self.view.camera.viewbox_mouse_event)
            self.view.camera._viewbox.events.mouse_release.connect(
                        self.view.camera.viewbox_mouse_event)
            self.view.camera._viewbox.events.mouse_wheel.connect(
                        self.view.camera.viewbox_mouse_event)



#================================= Event Functions Start ==================================#

    def on_key_press(self, event):
        # Set selection_flag and instruction text

        if event.text in self.selection_pool.keys():
            if not self.selection_flag:
                self.text.text = 'Now is %s selection mode, press %s to switch' % (self.selection_pool[event.text],
                                                                                   event.text)
                self.selection_flag = True
            else:
                self.text.text = 'Now is view mode, press %s to switch' % event.text
                self.selection_flag = False
            self.event_connect(self.selection_flag)
            self.selection_id = event.text
            # self.volume.visible = True

    def on_mouse_press(self, event):
        print('I wanna know mouse pos', event.pos)
        self.selection_origin = event.pos

        # Realize picking functionality and set origin mouse pos
        if event.button == 1 and self.selection_flag:
                # get the tr from canvas to the visual
                tr_back = self.volume.get_transform(map_from='canvas', map_to='visual')
                # tr_to_camera = self.view.camera.get_transform(map_from='canvas', map_to='camera')
                # camera_pos = tr_to_camera.map([400, 300])
                # print('really the camera pos?', camera_pos)
                cam = self.view.camera
                center_point = self.view.camera.center
                cam_point = self.view.camera.transform.map(center_point)

                trpos_start = np.insert(self.selection_origin, 2, 1)
                trpos_end = np.insert(self.selection_origin, 2, -1)
                print('before transpos start and end', trpos_start, trpos_end)

                trpos_start = tr_back.map(trpos_start)
                trpos_start = trpos_start[:3] / trpos_start[3]

                trpos_end = tr_back.map(trpos_end)
                trpos_end = trpos_end[:3] / trpos_end[3]
                print('after transpos start and end', trpos_start, trpos_end)

                pos = np.array([trpos_start, cam_point[:3]]) # CORRECT
                self.ray_line.set_data(pos=pos)
                self.rayline_text.text = 'start is %s \n end/camera is %s' % (str(np.round(trpos_start, 2)),
                                                                              str(np.round(cam_point[:3], 2)))
                self.canvas.update()

if __name__ == '__main__':
    appQt = QtGui.QApplication(sys.argv)
    view = DemoScene(keys='interactive')
    view.show()
    appQt.exec_()

__author__ = 'penny'

"""
    3D Ray Cursor based on Vispy library
    Date: May 21, 2016
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
        self.view.camera = scene.cameras.TurntableCamera(elevation=90,
                                                         azimuth=0,
                                                         fov=60,
                                                         center=(0, 0, 0),
                                                         distance=10.)

        # Data
        fitsdata = fits.open('test_data/l1448_13co.fits')
        self.vol_data = np.nan_to_num(fitsdata[0].data)
        self.vol_shape = self.vol_data.shape

        self.pos_data = np.indices(self.vol_shape).reshape(3, -1).transpose()

        self.volume = scene.visuals.Volume(self.vol_data,
                                           parent=self.view.scene)
        self.scale = (0.5, 0.5, 0.5)  # xyz
        self.trans = [-self.vol_shape[2]/2 * self.scale[0],
                      -self.vol_shape[1]/2 * self.scale[1],
                      -self.vol_shape[0]/2 * self.scale[2]]
        self.volume.transform = scene.STTransform(translate=self.trans,
                                                  scale=self.scale)

        # add some markers
        self.markers = scene.visuals.Markers(parent=self.view.scene)

        self.tr = self.volume.get_transform(map_from='visual', map_to='canvas')
        self.tr_back = self.volume.get_transform(map_from='canvas', map_to='visual')

        # Add a text instruction
        self.text = scene.visuals.Text('', color='red',
                                       pos=(self.canvas.size[0]/4.0,  20),
                                       parent=self.canvas.scene)

        # Max value pos text
        self.max_text = scene.visuals.Text('', color='yellow',
                                           parent=self.canvas.scene)

        # Add a cursor line
        self.ray_line = scene.visuals.Line(color='green', width=10,
                                           parent=self.view.scene)
        # Add text instruction for start point and end point
        self.rayline_text = scene.visuals.Text('', color='green', font_size=10,
                                               pos=(self.canvas.size[0]/2.0, 500),
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
        cam_events = self.view.camera._viewbox.events
        cam_mouse_event = self.view.camera.viewbox_mouse_event
        if flag:
            cam_events.mouse_move.disconnect(cam_mouse_event)
            cam_events.mouse_press.disconnect(cam_mouse_event)
            cam_events.mouse_release.disconnect(cam_mouse_event)
            cam_events.mouse_wheel.disconnect(cam_mouse_event)
        else:
            cam_events.mouse_move.connect(cam_mouse_event)
            cam_events.mouse_press.connect(cam_mouse_event)
            cam_events.mouse_release.connect(cam_mouse_event)
            cam_events.mouse_wheel.connect(cam_mouse_event)

    def on_key_press(self, event):
        # Set selection_flag and instruction text

        if event.text in self.selection_pool.keys():
            if not self.selection_flag:
                self.text.text = 'Now is %s selection mode, press %s to switch' \
                                 % (self.selection_pool[event.text], event.text)
                self.selection_flag = True
            else:
                self.text.text = 'Now is view mode, press %s to switch' \
                                 % event.text
                self.selection_flag = False
            self.event_connect(self.selection_flag)
            self.selection_id = event.text

    def on_mouse_press(self, event):
        """
        Use Markers to show the ray casting line, highlight max-value marker.
        """
        self.selection_origin = event.pos

        if event.button == 1 and self.selection_flag:
                pos = self.get_ray_line()
                
                # get intersected points in visual coordinate through z axis
                inter_pos = []

                for z in range(0, self.vol_shape[0]):
                    """
                      3D line defined with two points (x0, y0, z0) and (x1, y1, z1) as
                      (x - x1)/(x2 - x1) = (y - y1)/(y2 - y1) = (z - z1)/(z2 - z1) = t
                    """
                    z = z * self.scale[2] + self.trans[2]
                    t = (z - pos[0][2])/(pos[1][2] - pos[0][2])
                    x = t * (pos[1][0] - pos[0][0]) + pos[0][0]
                    y = t * (pos[1][1] - pos[0][1]) + pos[0][1]
                    inter_pos.append([x, y, z])

                inter_pos = np.array(inter_pos)

                # cut the line within the cube
                m1 = inter_pos[:, 0] > self.trans[0] # or =?  for x
                m2 = inter_pos[:, 0] < (self.vol_shape[2] * self.scale[0] + self.trans[0])
                m3 = inter_pos[:, 1] > self.trans[1]  # for y
                m4 = inter_pos[:, 1] < (self.vol_shape[1]*self.scale[1] + self.trans[1])
                inter_pos = inter_pos[m1 & m2 & m3 & m4]

                # set colors for markers
                colors = np.ones((inter_pos.shape[0], 4))
                colors[:] = (0.5, 0.5, 0, 1)

                inter_value = self.get_inter_value(inter_pos)  # value of intersected points
                assert inter_value.shape[0] == inter_pos.shape[0]

                # set max-value marker as different color
                max_data_value = inter_pos[np.argmax(inter_value), :3]
                if len(inter_value) != 0:
                    colors[np.argmax(inter_value)] = (1, 1, 0.5, 1)

                self.markers.set_data(pos=inter_pos, face_color=colors)

                # set text
                self.max_text.text = str(max_data_value)
                self.max_text.pos = self.selection_origin-5
                self.rayline_text.text = 'start is %s \n end/camera is %s' \
                                         % (str(np.round(inter_pos[0], 2)),
                                            str(np.round(inter_pos[-1], 2)))

                self.canvas.update()

    def get_inter_value(self, pos):
        """
        Get the value of intersected points.
        """
        inter_value = []
        for each_point in pos:
            x = (each_point[0] - self.trans[0])/self.scale[0]
            y = (each_point[1] - self.trans[1])/self.scale[1]
            z = (each_point[2] - self.trans[2])/self.scale[2]
            inter_value.append(self.vol_data[(z, y, x)])

        return np.array(inter_value)

    def get_ray_line(self):
        """
        Get the ray line from camera pos to the far point.

        :return: Start point and end point position.
        """
        tr_back = self.volume.get_transform(map_from='canvas', map_to='visual')

        center_point = self.view.camera.center
        cam_point = self.view.camera.transform.map(center_point)

        start = np.insert(self.selection_origin, 2, 1)
        start = tr_back.map(start)
        start = start[:3] / start[3]
        
        line_data = np.array([start, cam_point[:3]])
        # add a line to test the transform
        # test_ray_line = scene.visuals.Line(line_data, color='green',
                                           # width=10, parent=self.view.scene)
        self.view.add(test_ray_line)

        return np.array([start, cam_point[:3]])

if __name__ == '__main__':
    # print(self.canvas.scene.describe_tree(with_transform=True))
    appQt = QtGui.QApplication(sys.argv)
    view = DemoScene(keys='interactive')
    view.show()
    appQt.exec_()

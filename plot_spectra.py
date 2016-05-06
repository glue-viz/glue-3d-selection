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

        # Create canvas
        grid = self.canvas.central_widget.add_grid()

        # Second ViewBox uses a 3D perspective camera
        self.vb2 = grid.add_view(name='vb2', row=0, col=0, border_color='yellow')
        self.vb2.parent = self.canvas.scene
        self.vb2.camera = scene.TurntableCamera(elevation=90, azimuth=0, fov=60, center=(0,0,0), distance=10.)

        # Create two ViewBoxes, place side-by-side
        self.vb1 = grid.add_view(name='vb1', row=1, col=0, border_color='yellow')
        # First ViewBox uses a 2D pan/zoom camera
        self.vb1.camera = 'panzoom'
        self.spectra_line = scene.visuals.LinePlot(data=(0, 0), color='white', symbol='o', parent=self.vb1.scene)


        # add some axis
        # add some axes
        self.x_axis = scene.AxisWidget(orientation='bottom')
        # x_axis.stretch = (1, 0.1)
        grid.add_widget(self.x_axis, row=1, col=0)
        self.x_axis.link_view(self.vb1)
        self.y_axis = scene.AxisWidget(orientation='left')
        # y_axis.stretch = (0.1, 1)
        grid.add_widget(self.y_axis, row=1, col=0)
        self.y_axis.link_view(self.vb1)

        # Data
        fitsdata = fits.open('test_data/l1448_13co.fits')
        self.vol_data = np.nan_to_num(fitsdata[0].data)

        new_pos = np.transpose(self.vol_data)

        # TODO: replace the min&max threshold with real settings in Glue UI
        self.pos_data = np.indices(self.vol_data.shape).reshape(3,-1).transpose()

        self.volume = scene.visuals.Volume(self.vol_data, parent=self.vb2.scene)
        self.trans = [-self.vol_data.shape[2]/2, -self.vol_data.shape[1]/2, -self.vol_data.shape[0]/2]
        self.volume.transform = scene.STTransform(translate=self.trans)

        # add some markers
        self.markers = scene.visuals.Markers(parent=self.vb2.scene)

        self.tr = self.volume.get_transform(map_from='visual', map_to='canvas')
        self.tr_back = self.volume.get_transform(map_from='canvas', map_to='visual')

        # Add a text instruction
        self.text = scene.visuals.Text('', color='red', pos=(self.canvas.size[0]/4.0,  20), parent=self.canvas.scene)

        # Max value pos text
        self.max_text = scene.visuals.Text('', color='yellow', parent=self.canvas.scene)

        # Add a cursor line
        self.ray_line = scene.visuals.Line(color='green', width=10, parent=self.vb2.scene)
        # Add text instruction for start point and end point
        self.rayline_text = scene.visuals.Text('', color='green', font_size=10, pos=(self.canvas.size[0]/2.0, 500),
                                               parent=self.canvas.scene)

        # Add a 3D axis to keep us oriented
        axis = scene.visuals.XYZAxis(parent=self.vb2.scene)

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
            self.vb2.camera._viewbox.events.mouse_move.disconnect(
                    self.vb2.camera.viewbox_mouse_event)
            self.vb2.camera._viewbox.events.mouse_press.disconnect(
                        self.vb2.camera.viewbox_mouse_event)
            self.vb2.camera._viewbox.events.mouse_release.disconnect(
                        self.vb2.camera.viewbox_mouse_event)
            self.vb2.camera._viewbox.events.mouse_wheel.disconnect(
                        self.vb2.camera.viewbox_mouse_event)
        else:
            self.vb2.camera._viewbox.events.mouse_move.connect(
                    self.vb2.camera.viewbox_mouse_event)
            self.vb2.camera._viewbox.events.mouse_press.connect(
                        self.vb2.camera.viewbox_mouse_event)
            self.vb2.camera._viewbox.events.mouse_release.connect(
                        self.vb2.camera.viewbox_mouse_event)
            self.vb2.camera._viewbox.events.mouse_wheel.connect(
                        self.vb2.camera.viewbox_mouse_event)


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
                center_point = self.vb2.camera.center
                cam_point = self.vb2.camera.transform.map(center_point)

                trpos_start = np.insert(self.selection_origin, 2, 1)
                trpos_end = np.insert(self.selection_origin, 2, -1)
                print('before transpos start and end', trpos_start, trpos_end)

                trpos_start = tr_back.map(trpos_start)
                trpos_start = trpos_start[:3] / trpos_start[3]

                trpos_end = tr_back.map(trpos_end)
                trpos_end = trpos_end[:3] / trpos_end[3]
                print('after transpos start and end', trpos_start, trpos_end)

                pos = np.array([trpos_start, cam_point[:3]]) # CORRECT

                # finding the interected points position in visual coordinate through z axis
                inter_pos = np.zeros((self.vol_data.shape[0], 3))

                for z in range(-self.vol_data.shape[0]/2, self.vol_data.shape[0]/2):
                    #   3D line defined with two points (x0, y0, z0) and (x1, y1, z1) as
                    #   (x - x1)/(x2 - x1) = (y - y1)/(y2 - y1) = (z - z1)/(z2 - z1) = t
                    t = (z - pos[0][2])/(pos[1][2] - pos[0][2])
                    inter_pos[z, 0] = x = t * (pos[1][0] - pos[0][0]) + pos[0][0]
                    inter_pos[z, 1] = y = t * (pos[1][1] - pos[0][1]) + pos[0][1]
                    inter_pos[z, 2] = z

                # cut the line within the cube
                inter_pos = inter_pos[inter_pos[:, 0] > -self.vol_data.shape[2]/2]
                inter_pos = inter_pos[inter_pos[:, 0] < self.vol_data.shape[2]/2]
                inter_pos = inter_pos[inter_pos[:, 1] > -self.vol_data.shape[1]/2]
                inter_pos = inter_pos[inter_pos[:, 1] < self.vol_data.shape[1]/2]

                self.rayline_text.text = 'start is %s \n end/camera is %s' % (str(np.round(inter_pos[0], 2)),
                                                                              str(np.round(inter_pos[-1], 2)))
                # match the translate with volume data
                inter_pos = inter_pos.astype(int)-self.trans

                # set colors for markers visual
                colors = np.ones((inter_pos.shape[0], 4))
                colors[:] = (0.5, 0.5, 0, 1)

                inter_array = []
                for each_pos in inter_pos:
                    each_pos = tuple(each_pos[::-1])   # order of self.vol_data is [z, y, x]
                    inter_array.append(self.vol_data[each_pos])

                max_data_value = inter_pos[np.argmax(inter_array), :3]
                colors[np.argmax(inter_array)] = (1, 1, 0.5, 1)  # change the color of max value marker

                inter_pos = inter_pos.astype(int)+self.trans

                self.max_text.text = str(max_data_value)
                self.max_text.pos = self.selection_origin-5

                self.markers.set_data(pos=inter_pos, face_color=colors)

                # add the spectra plot
                x_lim = (min(inter_pos[:, 2])-self.trans[0], max(inter_pos[:, 2])-self.trans[0]) # using z range for x_lim
                y_lim = (min(inter_array), max(inter_array))

                self.x_axis.stretch = x_lim
                self.y_axis.stretch = y_lim

                line_pos = np.zeros((len(inter_array), 2), dtype=float)
                print('line_pos', line_pos)
                line_pos[:, 0] = inter_pos[:, 0]-self.trans[0]
                line_pos[:, 1] = inter_array

                # x axis is the z intervals, y axis is the point value
                self.spectra_line.set_data(line_pos)
                self.canvas.update()

if __name__ == '__main__':
    appQt = QtGui.QApplication(sys.argv)
    view = DemoScene(keys='interactive')
    view.show()
    appQt.exec_()

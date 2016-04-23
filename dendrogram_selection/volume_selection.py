__author__ = 'penny'

"""
    3D Points Picking & Selection tool based on Vispy library
    Date: Feb 15, 2016
    Usage: running the code with 'python selection_vispy_github.py'.
           press '1' or '2' to switch between view mode and selection mode for lasso method and picking method, respectively.
"""

import sys
import numpy as np

from PyQt4 import QtGui, QtCore
from vispy import app, scene
from matplotlib import path
import pyfits
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

from multivol import MultiVolume
from multivol import get_translucent_cmap

app.use_app('pyqt4')


def rectangle_vertice(center, height, width):
    # Borrow from _generate_vertices in vispy/visuals/rectangle.py

    half_height = height / 2.
    half_width = width / 2.

    bias1 = np.ones(4) * half_width
    bias2 = np.ones(4) * half_height

    corner1 = np.empty([1, 3], dtype=np.float32)
    corner2 = np.empty([1, 3], dtype=np.float32)
    corner3 = np.empty([1, 3], dtype=np.float32)
    corner4 = np.empty([1, 3], dtype=np.float32)

    corner1[:, 0] = center[0] - bias1[0]
    corner1[:, 1] = center[1] - bias2[0]
    corner1[:, 2] = 0


    corner2[:, 0] = center[0] + bias1[1]
    corner2[:, 1] = center[1] - bias2[1]
    corner2[:, 2] = 0


    corner3[:, 0] = center[0] + bias1[2]
    corner3[:, 1] = center[1] + bias2[2]
    corner3[:, 2] = 0


    corner4[:, 0] = center[0] - bias1[3]
    corner4[:, 1] = center[1] + bias2[3]
    corner4[:, 2] = 0

    # Get vertices between each corner of the rectangle for border drawing
    vertices = np.concatenate(([[center[0], center[1], 0.]],
                             [[center[0] - half_width, center[1], 0.]],
                             corner1,
                             [[center[0], center[1] - half_height, 0.]],
                             corner2,
                             [[center[0] + half_width, center[1], 0.]],
                             corner3,
                             [[center[0], center[1] + half_height, 0.]],
                             corner4,
                             [[center[0] - half_width, center[1], 0.]]))

    # vertices = np.array(output, dtype=np.float32)

    return vertices[1:, ..., :2]


def ellipse_vertice(center, radius, start_angle, span_angle, num_segments):
    # Borrow from _generate_vertices in vispy/visual/ellipse.py

    if isinstance(radius, (list, tuple)):
        if len(radius) == 2:
            xr, yr = radius
        else:
            raise ValueError("radius must be float or 2 value tuple/list"
                             " (got %s of length %d)" % (type(radius),
                                                         len(radius)))
    else:
        xr = yr = radius

    start_angle = np.deg2rad(start_angle)

    vertices = np.empty([num_segments + 2, 2], dtype=np.float32) # Segment as 1000

    # split the total angle into num_segments intances
    theta = np.linspace(start_angle,
                        start_angle + np.deg2rad(span_angle),
                        num_segments + 1)

    # PolarProjection
    vertices[:-1, 0] = center[0] + xr * np.cos(theta)
    vertices[:-1, 1] = center[1] + yr * np.sin(theta)

    # close the curve
    vertices[num_segments + 1] = np.float32([center[0], center[1]])

    return vertices[:-1]


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
        self.canvas.events.mouse_release.connect(self.on_mouse_release)
        self.canvas.events.mouse_move.connect(self.on_mouse_move)
        self.canvas.events.key_press.connect(self.on_key_press)

        # Setup some defaults
        self.mesh = None
        self.selected = []
        self.white = (1.0, 1.0, 1.0, 1.0)
        self.black = (0.0, 0.0, 0.0, 0.0)

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

        self.tr = self.volume.node_transform(self.view)  # ChainTransform

        # Add a text instruction
        self.text = scene.visuals.Text('', color='white', pos=(self.canvas.size[0]/4.0,  20), parent=self.canvas.scene)

        # Add a 3D axis to keep us oriented
        axis = scene.visuals.XYZAxis(parent=self.view.scene)

        # Set up for lasso drawing
        self.line_pos = []
        self.line = scene.visuals.Line(color='yellow', method='gl', parent=self.canvas.scene)

        # Selection
        self.selection_flag = False
        self.selection_pool = {'1': 'lasso', '2': 'rectangle', '3': 'ellipse', '4': 'pick'}
        self.selection_id = '1'  # default as 1
        self.selection_origin = (0, 0)

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

    def mark_selected(self):
        # Change the color of the selected point
        reds = get_translucent_cmap(1, 0, 0)

        select_data = np.transpose(self.vol_data)
        not_select = np.logical_not(self.selected)
        np.place(select_data, not_select, 0)
        select_data = np.transpose(select_data)
        
        print('select_data is', select_data, select_data.shape)
        maxpos = np.unravel_index(select_data.argmax(), select_data.shape)
        print('got the max pos', maxpos)
        self.volume_pool.append((select_data, (1, 6), reds))

        # TODO: no set_data function available in multi_volume_visual
        self.volume._update_all_volumes(self.volume_pool)
        print('self.volume_pool', len(self.volume_pool))
        self.canvas.update()

    def on_key_press(self, event):
        # Set selection_flag and instruction text

        if event.text in self.selection_pool.keys():
            if not self.selection_flag:
                self.text.text = 'Now is %s selection mode, press %s to switch' % (self.selection_pool[event.text], event.text)
                self.selection_flag = True
            else:
                self.text.text = 'Now is view mode, press %s to switch' % event.text
                self.selection_flag = False
            self.event_connect(self.selection_flag)
            self.selection_id = event.text
            self.volume.visible = True

    def on_mouse_press(self, event):
        # Realize picking functionality and set origin mouse pos

        if event.button == 1 and self.selection_flag:
            if self.selection_id == '4':
                # Ray intersection on the CPU to highlight the selected point(s)
                data = self.tr.map(self.pos_data)[:, :2]  # Map coordinates
                print('data after tr.map', data)
                m1 = data > (event.pos - 4)
                m2 = data < (event.pos + 4)

                pick_selected = np.argwhere(m1[:,0] & m1[:,1] & m2[:,0] & m2[:,1])
                len_mask = self.vol_data.shape[0]*self.vol_data.shape[1]*self.vol_data.shape[2]
                
                full_mask = np.zeros(len_mask)
                full_mask[pick_selected] = True
                self.selected = full_mask
                print('self.selected is', self.selected, len(self.selected))
                self.mark_selected()

            else:
                self.selection_origin = event.pos

    def on_mouse_release(self, event):
        # Identify selected points and mark them

        if event.button == 1 and self.selection_flag and self.selection_id is not '4':
            data = self.tr.map(self.pos_data)[:, :2]

            if self.selection_id in ['1', '2', '3']:
                selection_path = path.Path(self.line_pos, closed=True)
                mask = selection_path.contains_points(data)

                self.selected = mask
                print('mask len', len(mask), mask)
                self.mark_selected()

                # Reset lasso
                self.line_pos = []  # TODO: Empty pos input is not allowed for line_visual
                self.line.set_data(np.array(self.line_pos))
                self.line.update()

            if self.selection_id in ['2', '3']:
                self.selection_origin = None

    def on_mouse_move(self, event):
        # Draw lasso/rectangle/ellipse shape with mouse dragging

        if event.button == 1 and event.is_dragging and self.selection_flag:
            if self.selection_id == '1':
                self.line_pos.append(event.pos)
                self.line.set_data(np.array(self.line_pos))

            if self.selection_id in ['2', '3']:
                width = event.pos[0] - self.selection_origin[0]
                height = event.pos[1] - self.selection_origin[1]
                center = (width/2. + self.selection_origin[0], height/2.+self.selection_origin[1], 0)

                if self.selection_id == '2':
                    self.line_pos = rectangle_vertice(center, height, width)
                    self.line.set_data(np.array(self.line_pos))

                if self.selection_id == '3':
                    self.line_pos = ellipse_vertice(center, radius=(np.abs(width/2.), np.abs(height/2.)),
                                                    start_angle=0., span_angle=360., num_segments=500)
                    self.line.set_data(pos=np.array(self.line_pos), connect='strip')


if __name__ == '__main__':
    appQt = QtGui.QApplication(sys.argv)
    view = DemoScene(keys='interactive')
    view.show()
    appQt.exec_()

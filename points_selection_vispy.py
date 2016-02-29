"""
    3D Points Picking & Selection tool based on Vispy library
    Author: Penny Xuran Qian (NAOC/SAO) pennyqxr@gmail.com
    Date: Feb 15, 2016
    License: GPL

    Usage: running the code with 'python selection_vispy_github.py'.
           press '1' or '2' to switch between view mode and selection mode for lasso method and picking method, respectively.
"""

import sys
import numpy as np

from PyQt4 import QtGui, QtCore
from vispy import app, scene
from matplotlib import path

app.use_app('pyqt4')


def rectangle_vertice(center, height, width):
    # Borrow from _generate_vertices in vispy/visual/rectangle.py

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

    output = np.concatenate(([[center[0], center[1], 0.]],
                             [[center[0] - half_width, center[1], 0.]],
                             corner1,
                             [[center[0], center[1] - half_height, 0.]],
                             corner2,
                             [[center[0] + half_width, center[1], 0.]],
                             corner3,
                             [[center[0], center[1] + half_height, 0.]],
                             corner4,
                             [[center[0] - half_width, center[1], 0.]]))

    vertices = np.array(output, dtype=np.float32)
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
        self.resize(500,500)
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

        # Test data
        self.data = np.random.uniform(-1, 1, size=(200, 3))
        self.minz = np.min(self.data[:,2])
        self.facecolor = np.ones((200,4), dtype=np.float)
        self.ptsize = 10
        self.scatter = scene.visuals.Markers()
        self.scatter.set_data(self.data, face_color=self.facecolor, edge_color=None, size=self.ptsize)
        self.view.add(self.scatter)

        # Add a text instruction
        self.text = scene.visuals.Text('', color='white', pos=(self.canvas.size[0]/4.0,  20), parent=self.canvas.scene)

        self.tr = self.scatter.node_transform(self.view)

        # Add a 3D axis to keep us oriented
        axis = scene.visuals.XYZAxis(parent=self.view.scene)

        # Set up for lasso drawing
        self.line_pos = []
        self.line = scene.visuals.Line(color='white', method='gl', parent=self.canvas.scene)

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

        self.facecolor[self.facecolor[:,1] != 1.0] = self.white
        self.scatter.set_data(self.data, face_color=self.facecolor, size=self.ptsize)
        for i in self.selected:
            self.facecolor[i] = [1.0, 0.0, 0.0, 1]

        self.scatter.set_data(self.data, face_color = self.facecolor, size=self.ptsize)
        self.scatter.update()

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

    def on_mouse_press(self, event):
        # Realize picking functionality and set origin mouse pos

        if event.button == 1 and self.selection_flag:
            if self.selection_id == '4':
                # Ray intersection on the CPU to highlight the selected point(s)

                data = self.tr.map(self.data)[:, :2]
                m1 = data > (event.pos - 4)
                m2 = data < (event.pos + 4)

                self.selected = np.argwhere(m1[:,0] & m1[:,1] & m2[:,0] & m2[:,1])
                self.mark_selected()

            else:
                self.selection_origin = event.pos

    def on_mouse_release(self, event):
        # Identify selected points and mark them

        if event.button == 1 and self.selection_flag and self.selection_id is not '4':
            self.facecolor[self.facecolor[:,1] != 1.0] = self.white

            data = self.tr.map(self.data)[:,:2]

            if self.selection_id in ['1', '2', '3']:
                selection_path = path.Path(self.line_pos, closed=True)
                mask = [selection_path.contains_points(data)]

                self.selected = mask
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

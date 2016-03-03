__author__ = 'penny'

import numpy as np
from vispy import app, scene
from matplotlib import path


class SelectionCommon(object):
    def __init__(self, canvas, view, data, scatter, facecolor):
        self._canvas = canvas
        self._view = view

        # Initialize drawing visual
        self.line_pos = []
        self.line = scene.visuals.Line(color='white', method='gl', parent=self._canvas.scene)

        # Selection defaults
        self.selection_flag = False
        self.selection_pool = {'1': 'lasso', '2': 'rectangle', '3': 'ellipse', '4': 'pick'}
        self.selection_id = '1'  # default as 1
        self.selection_origin = (0, 0)
        self.selected = []

        self.white = (1.0, 1.0, 1.0, 1.0)
        self.black = (0.0, 0.0, 0.0, 0.0)
        self.facecolor = facecolor

        # Scatter plot data, visual and projection
        self._data = data
        self._scatter = scatter
        self.tr = self._scatter.node_transform(self._view)

        # Add a text instruction
        self.text = scene.visuals.Text('', color='white', pos=(self._canvas.size[0] / 4.0, 20),
                                       parent=self._canvas.scene)

        self.canvas_event_connect()

    def canvas_event_connect(self):
        self._canvas.events.mouse_press.connect(self.on_mouse_press)
        self._canvas.events.mouse_release.connect(self.on_mouse_release)
        self._canvas.events.mouse_move.connect(self.on_mouse_move)
        self._canvas.events.key_press.connect(self.on_key_press)

    def camera_viewbox_events(self, flag):
        """
        Disconnect camera events with viewbox when drawing selection line
        """
        if flag:
            self._view.camera._viewbox.events.mouse_move.disconnect(
                self._view.camera.viewbox_mouse_event)
            self._view.camera._viewbox.events.mouse_press.disconnect(
                self._view.camera.viewbox_mouse_event)
            self._view.camera._viewbox.events.mouse_release.disconnect(
                self._view.camera.viewbox_mouse_event)
            self._view.camera._viewbox.events.mouse_wheel.disconnect(
                self._view.camera.viewbox_mouse_event)
        else:
            self._view.camera._viewbox.events.mouse_move.connect(
                self._view.camera.viewbox_mouse_event)
            self._view.camera._viewbox.events.mouse_press.connect(
                self._view.camera.viewbox_mouse_event)
            self._view.camera._viewbox.events.mouse_release.connect(
                self._view.camera.viewbox_mouse_event)
            self._view.camera._viewbox.events.mouse_wheel.connect(
                self._view.camera.viewbox_mouse_event)

    def rectangle_vertice(self, center, height, width):
        """
        Borrow from _generate_vertices in vispy/visuals/rectangle.py
        """
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

    def ellipse_vertice(self, center, radius, start_angle, span_angle, num_segments):
        """
        Borrow from _generate_vertices in vispy/visual/ellipse.py
        """
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

        vertices = np.empty([num_segments + 2, 2], dtype=np.float32)  # Segment as 1000

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

    def mark_selected(self):
        self.facecolor[self.facecolor[:, 1] != 1.0] = self.white
        self._scatter.set_data(self._data, face_color=self.facecolor)
        for i in self.selected:
            self.facecolor[i] = [1.0, 0.0, 0.0, 1]

        self._scatter.set_data(self._data, face_color=self.facecolor)
        self._scatter.update()

    # TODO: replace with toolbar handle
    def on_key_press(self, event):
        """
        Set selection_flag and instruction text
        """
        if event.text in self.selection_pool.keys():
            if not self.selection_flag:
                self.text.text = 'Now is %s selection mode, press %s to switch' % (
                    self.selection_pool[event.text], event.text)
                self.selection_flag = True
            else:
                self.text.text = 'Now is view mode, press %s to switch' % event.text
                self.selection_flag = False
            self.camera_viewbox_events(self.selection_flag)
            self.selection_id = event.text

    def on_mouse_press(self, event):
        """
        Realize picking functionality and set origin mouse pos
        """
        if event.button == 1 and self.selection_flag:
            if self.selection_id == '4':

                # Ray intersection on the CPU to highlight the selected point(s)
                data = self.tr.map(self._data)[:, :2]
                m1 = data > (event.pos - 4)
                m2 = data < (event.pos + 4)

                self.selected = np.argwhere(m1[:, 0] & m1[:, 1] & m2[:, 0] & m2[:, 1])
                self.mark_selected()

            else:
                self.selection_origin = event.pos

    def on_mouse_release(self, event):
        """
        Identify selected points and mark them
        """
        if event.button == 1 and self.selection_flag and self.selection_id is not '4':
            self.facecolor[self.facecolor[:, 1] != 1.0] = self.white

            data = self.tr.map(self._data)[:, :2]

            if self.selection_id in ['1', '2', '3']:
                selection_path = path.Path(self.line_pos, closed=True)
                mask = [selection_path.contains_points(data)]

                self.selected = mask
                self.mark_selected()

                # Reset selection line
                self.line_pos = []  # TODO: Empty pos input is not allowed for line_visual
                self.line.set_data(np.array(self.line_pos))
                self.line.update()

            if self.selection_id in ['2', '3']:
                self.selection_origin = None

    def on_mouse_move(self, event):
        """
        Draw selection line along dragging mouse
        """
        if event.button == 1 and event.is_dragging and self.selection_flag:
            if self.selection_id == '1':
                self.line_pos.append(event.pos)
                self.line.set_data(np.array(self.line_pos))

            if self.selection_id in ['2', '3']:
                width = event.pos[0] - self.selection_origin[0]
                height = event.pos[1] - self.selection_origin[1]
                center = (width / 2. + self.selection_origin[0], height / 2. + self.selection_origin[1], 0)

                if self.selection_id == '2':
                    self.line_pos = self.rectangle_vertice(center, height, width)
                    self.line.set_data(np.array(self.line_pos))

                if self.selection_id == '3':
                    self.line_pos = self.ellipse_vertice(center, radius=(np.abs(width / 2.), np.abs(height / 2.)),
                                                         start_angle=0., span_angle=360., num_segments=500)
                    self.line.set_data(pos=np.array(self.line_pos), connect='strip')
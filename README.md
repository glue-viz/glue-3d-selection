Experimental 3D selection scripts 
================================

About
------

For scatter_selection folder, there are three files as: `points_selection_vispy.py`(a standalone script), `selection_common.py` and `selection_widget.py`(these two are functionally diveded from `points_selection_vispy.py`).

For volume_selection folder, it's mostly the same as scatter points case, except using `MultiVolume` developed by @astrofrog on the [repository](https://github.com/astrofrog/vispy-multivol).

How to use 
-----------
Users can either run the standalone script with
`python points_selection_vispy.py` or `python volume_selection.py`

or run the widget scripy with
`python selection_widget.py`

There are four selection modes switched with key press:
* '1' stands for free lasso selection;
* '2' stands for rectangle shape selection;
* '3' stands for ellipse shape selection;
* '4' stands for picking (not available for `volume selection` yet)

While under selection mode, the mouse operation for controlling the viewport is not available, users could press one of above key press to swith between selection & view mode.


from pyqtgraph.Qt import QtGui
import pyqtgraph as pg
import numpy as np


class Scope:
    def __init__(self, title, min_y, max_y):
        self.refresh_counter = 0
        self.win = pg.GraphicsWindow(title=title) # creates a window
        p = self.win.addPlot()  # creates empty space for the plot in the window
        self.curve1 = p.plot(pen=pg.mkPen('g', width=2))    # create an empty "plot" (a curve to plot)
        self.curve2 = p.plot(pen=pg.mkPen('r', width=2))    # create an empty "plot" (a curve to plot)
        n_samples = 1000                       # width of the window displaying the curve
        self.Xm = np.linspace(0, 0, n_samples)          # create array that will contain the relevant time series
        self.Xm2 = np.linspace(0, 0, n_samples)          # create array that will contain the relevant time series
        self.time = np.linspace(0, 0, n_samples)
        self.vb1 = self.curve1.getViewBox()
        self.vb2 = self.curve2.getViewBox()
        self.vb1.enableAutoRange(axis='y', enable=False)
        self.vb2.enableAutoRange(axis='y', enable=False)
        p.setYRange(min_y, max_y)

    def update(self, time, y1, y2):
        self.Xm[:-1] = self.Xm[1:]                      # shift data in the temporal mean 1 sample left
        self.Xm[-1] = y1                 # vector containing the instantaneous values
        self.Xm2[:-1] = self.Xm2[1:]                      # shift data in the temporal mean 1 sample left
        self.Xm2[-1] = y2                 # vector containing the instantaneous values
        self.time[:-1] = self.time[1:]
        self.time[-1] = time
        self.curve1.setData(self.time, self.Xm)                     # set the curve with this data
        self.curve2.setData(self.time, self.Xm2)                     # set the curve with this data
        self.refresh_counter += 1
        if self.refresh_counter == 100:
            self.refresh_counter = 0
            QtGui.QApplication.processEvents()

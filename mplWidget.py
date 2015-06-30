from PyQt4 import QtGui
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure

class MplCanvas(FigureCanvas):
	def __init__(self):
		self.fig = Figure()
		self.ax  = self.fig.add_subplot(111)
		FigureCanvas.__init__(self, self.fig)

class matplotlibWidget(QtGui.QWidget):
	def __init__(self,parent=None):
		QtGui.QWidget.__init__(self,parent)
		self.canvas = MplCanvas()
		self.toolbar = NavigationToolbar(self.canvas,self)
		self.vbl    = QtGui.QVBoxLayout()
		self.vbl.addWidget(self.canvas)
		self.vbl.addWidget(self.toolbar)
		self.setLayout(self.vbl)

	
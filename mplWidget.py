from PyQt4 import QtGui
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib import __version__
from matplotlib.figure import Figure
truncv = __version__[0:3]
if truncv=='1.5':
	from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
else: # earlier versions of Matplotlib OK
	from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar


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

	


try:
    from PyQt4.QtGui import QWidget, QVBoxLayout
   
except ImportError:
    try:
        print("PyQt4 not found. Using PyQt5.")
        from PyQt5.QtWidgets import QWidget
        from PyQt5.QtWidgets import QVBoxLayout
    except ImportError:
        raise Exception("PyQt5 not found. Please install either PyQt4 or PyQt5.")

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib import __version__

truncv = __version__[0:3]
if truncv=='1.5':
	from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
else: # earlier versions of Matplotlib OK
	try:
		from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
	except ImportError:
		from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

class MplCanvas(FigureCanvas):
	def __init__(self):
		self.fig = Figure()
		self.ax  = self.fig.add_subplot(111)
		FigureCanvas.__init__(self, self.fig)

class matplotlibWidget(QWidget):
	def __init__(self,parent=None):
		QWidget.__init__(self,parent)
		self.canvas = MplCanvas()
		self.toolbar= NavigationToolbar(self.canvas,self)
		self.vbl    = QVBoxLayout()
		self.vbl.addWidget(self.canvas)
		self.vbl.addWidget(self.toolbar)
		self.setLayout(self.vbl)

	

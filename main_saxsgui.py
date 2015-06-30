
# 06/10/2015 - added back Fortran module for making smeared transformation matrix
#            - the P(r) saved from the window is integrated to 1 for easier overlay later
# 07/16/2014 - added Error Bar display support for log-scale, since pyqtgraph doesn't have this
#              this was done by sub-classing ErrorBarItem (see below)
# 07/15/2014 - took out automatic Guinier analysis, code is buggy and needs to be fixed

from PyQt4.QtCore import *
from PyQt4.QtGui import *
import sys
import os
import time
from saxsgui import Ui_SAXSgui
import pyqtgraph as pg
from saxsmod import saxsdata, beamprofile, fitline, iftv2, grideval
from scipy.integrate import simps
from numpy import exp, log, log10, array, loadtxt, linspace, zeros,\
                 unravel_index,pi, dot, sqrt, arange, sin
from scipy.optimize import minimize

class TextItem(pg.TextItem):

    def __init__(self, text, **kwargs):
        pg.TextItem.__init__(self, text, **kwargs)
        self.logMode = [False, False]
        self.savex = self.pos()[0]
        self.savey = self.pos()[1]

    def setLogMode(self, xMode, yMode):
        if self.logMode == [xMode, yMode]:
            return
        self.logMode = [xMode, yMode]
        self.update_pos()

    def update_pos(self):
        x = self.savex
        y = self.savey
        if self.logMode[0]:
            x = log10(x)
        if self.logMode[1]:
            y = log10(y)
        self.setPos(x,y)


# ErrorBar support for ErrorBarItem (pyqtgraph doesn't currently have one yet)
class ErrorBarItem(pg.ErrorBarItem):

    def __init__(self, **opts):
        pg.ErrorBarItem.__init__(self, **opts)
        self.opts['logMode'] = [False, False]
        # keep original values
        self.saveopts = dict(x=self.opts['x'],\
                             y=self.opts['y'],\
                             width=self.opts['width'],\
                             height=self.opts['height'])

    def setLogMode(self, xMode, yMode):
        if self.opts['logMode'] == [xMode, yMode]:
            return
        self.opts['logMode'] = [xMode, yMode]
        self.updateData()
        self.informViewBoundsChanged()

    def updateData(self, **opts):
        # this dict update handles new data range being displayed
        self.opts.update(**opts)
        self.saveopts.update(**opts)

        x = self.saveopts['x']
        y = self.saveopts['y']
        h = self.saveopts['height']
        w = self.saveopts['width']

        if self.opts['logMode'][0]: # x log mode
            if w is not None:
                self.opts['width'] = w/x
            self.opts['x'] = log10(x)

        else:
            self.opts['x'] = x
            self.opts['width']= w

        if self.opts['logMode'][1]: # y log mode
            if h is not None: 
                self.opts['height'] = h/y
            self.opts['y'] = log10(y)

        else:
            self.opts['y'] = y
            self.opts['height'] = h

        self.drawPath() 

pg.setConfigOption('background','w')
pg.setConfigOption('foreground','k')
pg.setConfigOption('antialias',True)

class saxsgui_mainwindow(Ui_SAXSgui):
    def __init__(self):
        self.mainWindow = QMainWindow()
        self.setupUi(self.mainWindow)
        self.mainWindow.setWindowTitle("UCSF SAXS v0.5")
        # set connections to GUI here

        # initialize GUI with plots
        self.initializeGUI()
        self.initializeparams()

        # set up shortcuts for keyboard
        self.actionOpen.setShortcut(QKeySequence.Open)
        self.actionClose.setShortcut(QKeySequence.Quit)

        self.mainWindow.connect(self.actionOpen, SIGNAL("triggered()"), self.openfile)
        self.mainWindow.connect(self.checkLogY, SIGNAL("stateChanged(int)"), self.axisupdate)
        self.mainWindow.connect(self.checkLogX, SIGNAL("stateChanged(int)"), self.axisupdate)
        self.mainWindow.connect(self.actionLoad_Beam_Profile, SIGNAL("triggered()"),self.openbeamfile)
        self.mainWindow.connect(self.solvepr_button, SIGNAL("clicked()"),self.runift)
        self.mainWindow.connect(self.calcbiftgrid, SIGNAL("clicked()"),self.biftgrideval)
        self.mainWindow.connect(self.refinebift, SIGNAL("clicked()"),self.biftsimplex)
        self.mainWindow.connect(self.runprimaryanalysis, SIGNAL("clicked()"),self.manualGuinier)
        self.mainWindow.connect(self.nskip_spinbox, SIGNAL("valueChanged(int)"),self.adjustrange)
        self.mainWindow.connect(self.nskip2_spinbox, SIGNAL("valueChanged(int)"),self.adjustrange)
        self.mainWindow.connect(self.actionClose, SIGNAL("triggered()"), self.mainWindow.close)
        self.mainWindow.connect(self.actionSave_as_GNOM_file, SIGNAL("triggered()"),self.saveGNOMfile)

    def initializeparams(self):
        self.filename = None
        self.beamname = None
        self.cwd = os.getenv("HOME")
        self.supportedFormats = ("dat","pdh","txt","the")
        self.data = None
        self.statusbar.showMessage("Ready")     

    def initializeGUI(self):
        x = arange(0.0008, 3.14*8, 0.002)
        y = sin(x)/x

        # set form validation on input QLineEdits
        float_validator = QDoubleValidator()
        int_validator   = QIntValidator()
        self.dmax_input.setValidator(float_validator)
        self.logalpha_input.setValidator(float_validator)
        self.Nr_input.setValidator(int_validator)

        self.dmax_input.setText(str(80.0))
        self.logalpha_input.setText(str(6.0))
        self.Nr_input.setText(str(120))

        # set up & initialize plotting of raw data 
        self.redpen = pg.mkPen('#d46a6a',width=2)
        self.blackpen = pg.mkPen('k',width=2)
        self.bluepen = pg.mkPen('#6495ed',width=2)
        self.graydashedpen = pg.mkPen('#a9a9a9',width=2,style=Qt.DashLine)
        
        self.rawplot = pg.PlotItem(title="Raw Data",name="rawplot")
        self.rawplot.setLabel("bottom","Scattering Angle (q)")
        self.rawplot.setLabel("left","I(q) (arbitrary)")
        self.rawplot.enableAutoRange()
        self.fourierView.addItem(self.rawplot)
        self.rawplot.plot(x,y,pen=self.redpen)

        # Real Space Plot
        self.prplot = pg.PlotItem(title="P(r) solution",name="prplot")
        self.prplot.setLabel("bottom","Interatomic Distance")
        self.prplot.setLabel("left","P(r) (arbitrary units)")
        self.prplot.enableAutoRange()
        self.realspaceView.addItem(self.prplot)

        # Beam Profile plot
        self.beamplot = pg.PlotItem(title="Beam Profile",name="beamplot")
        self.beamplot.enableAutoRange()
        self.beamview.addItem(self.beamplot)

        self.guinierplot = pg.PlotItem(title="Guinier Analysis", name="GuinierPlot")
        self.kratkyplot  = pg.PlotItem(title="Normalized Kratky Plot", name="KratkyPlot")
        self.primaryanalysis_view.addItem(self.guinierplot,row=1,col=1)
        self.primaryanalysis_view.addItem(self.kratkyplot,row=1,col=2)
        
    def axisupdate(self):
        ylog = self.checkLogY.checkState() > 0
        xlog = self.checkLogX.checkState() > 0
        
        self.rawplot.setLogMode(x=xlog,y=ylog)

    def openfile(self):
        cwd = (os.path.dirname(self.filename)\
                if self.filename is not None else self.cwd)
        formats = (["*.{0}".format(unicode(formats).lower())\
                for formats in self.supportedFormats])
        fname = unicode(QFileDialog.getOpenFileName(self.mainWindow,\
                                                    "Open raw SAXS data ...",\
                                                    cwd,
                                                    "SAXS Data ({0})".\
                                                    format(" ".join(formats))))
        if fname:
            self.filename = fname
            self.loadfile(fname)

    def openbeamfile(self):
        cwd = (os.path.dirname(self.beamname)\
                if self.beamname is not None else self.cwd)
        formats = (["*.{0}".format(unicode(formats).lower())\
                for formats in self.supportedFormats])
        fname = unicode(QFileDialog.getOpenFileName(self.mainWindow,\
                                                    "Open Beam Profile ...",\
                                                    cwd,
                                                    "Beam Profile ({0})".\
                                                    format(" ".join(formats))))
        if fname:
            self.beamname = fname
            beam_fname = fname.split('/')[-1]
            self.beamdata = beamprofile(self.beamname)
            self.beamdata.normalize()
            self.y = self.beamdata.y
            self.Wy= self.beamdata.Wy
            self.beamplot.clear()
            self.beamplot.setLabels(title=beam_fname)
            self.beamplot.plot(self.y, self.Wy, pen=self.redpen)

    def saveGNOMfile(self):

        cwd = (os.path.dirname(self.filename)\
                if self.filename is not None else self.cwd)

        fname = unicode(QFileDialog.getSaveFileName(self.mainWindow,\
                                                    "Save as GNOM file ...",\
                                                    cwd))
        self.data.writeGNOM(fname)

    def loadfile(self,filename):
        self.data = saxsdata(self.filename)       
        raw_fname = self.filename.split('/')[-1]
        self.rawplot.clear()

        q  = self.data.q
        Iq = self.data.Iq
        sd = self.data.sd
        self.data.solved = False

        qmin = q.min()
        qmax = q.max()

        # set Qmin/max label to inform user angular range
        self.qmin_label.setText("{0:8.6f},(pi/q_min {1:6.2f})".format(qmin,pi/qmin))
        self.qmax_label.setText("{0:8.6f}".format(qmax))
        # set the guinier regime spinbox (upper limit)
        self.ub_guinier_spinbox.setMaximum(len(q))

        self.data.datarange = [0, len(q)]
        
        if qmax<1:
            increment = 0.001
        if qmax>1:
            increment = 0.01

        self.datarange = [0,len(q)]
        # reset the QSpinBoxes
        self.nskip_spinbox.setValue(0)
        self.nskip2_spinbox.setValue(0)

        # set errorbar for each point
        self.Iq_errorbar = ErrorBarItem(x=q,y=Iq,height=sd,beam=0,pen={'color':'#a9a9a9','width':2})

        # plot data 
        self.rawplot.clear()
        self.rawplot.setLabels(title=raw_fname)
        self.rawplot.addItem(self.Iq_errorbar)
        self.rawplot.plot(q,Iq,pen=self.redpen)


    def adjustrange(self):

        qmin_trunc = int(self.nskip_spinbox.value())
        qmax_trunc = len(self.data.q) - int(self.nskip2_spinbox.value())
        try :
            self.data.datarange = [qmin_trunc, qmax_trunc]
            wrkq = self.data.q[qmin_trunc:qmax_trunc]
            self.qmin_label.setText("{0:8.6f},(pi/q_min {1:6.2f})".format(wrkq.min(),pi/wrkq.min()))
            self.qmax_label.setText("{0:8.6f}".format(wrkq.max()))
            self.updateplot()
        except:
            raise

    def updateplot(self):

        qmin_id = self.data.datarange[0]
        qmax_id = self.data.datarange[1]

        self.rawplot.clear()
        self.Iq_errorbar.updateData(x=self.data.q[qmin_id:qmax_id],\
                                 y=self.data.Iq[qmin_id:qmax_id],\
                                 height=self.data.sd[qmin_id:qmax_id])
        self.rawplot.addItem(self.Iq_errorbar)
        self.rawplot.plot(self.data.q[qmin_id:qmax_id],\
                          self.data.Iq[qmin_id:qmax_id],\
                          pen=self.redpen)


    def getiftparams(self):
        self.data.Dmax = float(self.dmax_input.text())
        self.data.logalpha = float(self.logalpha_input.text())
        self.data.Nr     = int(self.Nr_input.text())

        weighdata = self.varianceweighted.checkState()
        smeared   = self.smeared.checkState()

        if weighdata==0:
            self.weighdata = False
        elif weighdata>0:
            self.weighdata = True

        if smeared==0:
            self.data.smeared = False
        elif weighdata>0:
            self.data.smeared = True


    def runift(self):
        # obtain parameters for solving P(r)
        self.getiftparams()

        # run fortran wrapper ift from iftreg module
        alpha = exp(self.data.logalpha)
        Dmax  = self.data.Dmax
        Nr = self.data.Nr

        if self.data.smeared:
            y  = self.y
            Wy = self.Wy
        else:
            y = None
            Wy = None
            
        q  = self.data.q
        Iq = self.data.Iq
        sd = self.data.sd
        
        qmin_id = self.data.datarange[0]
        qmax_id = self.data.datarange[1]

        q  = q[qmin_id:qmax_id]
        Iq = Iq[qmin_id:qmax_id]
        sd = sd[qmin_id:qmax_id]

        Jreg,Ireg, Jreg_extrap, Ireg_extrap, q_full, r,pr,evi = iftv2(alpha,Dmax,q,Iq,sd,Nr,y,Wy,self.weighdata,self.data.smeared)

        chi = ((Iq-Jreg)**2/sd**2).mean()

        # assign solution to saxsdata
        self.data.r, self.data.pr = r, pr
        self.data.Jreg = Jreg
        self.data.Ireg = Ireg
        self.data.solved = True

        # calculate real-space Rg, etc.
        dr = r[1]-r[0]
        Rg_sq = 0.5 * (pr * r**2 * dr).sum() / (pr*dr).sum()
        self.data.Rg = sqrt(Rg_sq)

        # Calculate SAXSMoW ...
        dq = q[1]-q[0]
        Iq_q = Ireg_extrap * q_full 
        # Integrate using Trapezoid rule
        #Qinv = 0.5*dq * (Iq_q[0] + Iq_q[-1] + 2*(Iq_q[1:-1]).sum())
        # Integrate using Simpsons's rule
        Qinv  = simps(y=Iq_q,x=q_full)

        Vc = Ireg_extrap[0]/Qinv
        
        ramboratio = Vc**2 / self.data.Rg

        # from saxsmow calibration
        # slope, offset  = saxsmow(q.max())
        
        # using protein density 1.35 g/mL
        # 1.35 g/cm^3 x 10^-3 kg/g x 1 Dalton/1.66e-27kg x 1cm^3/1e24 Angstrom^3
        # is 0.81325 Dalton/Angstrom^3
        #saxsmw = (slope*appVol + offset) * 0.81325e-3 
        # Using the numbers from biosis.net (Robert Rambo)
        # ONLY for protein / complexes, will have to be different for nucleic acid / complex
        saxsmw = (ramboratio/0.1231)**1.0
        saxsmw = saxsmw * 0.001
        chi_str = "Chi^2: {0:5.3f}\nI(0): {1:8.3E}\nM.W.: {2:4.1f}kDa".format(chi,Ireg_extrap[0],saxsmw)
        chi_textitem = TextItem(chi_str, color='k')
        chi_textitem.setPos(q.max()*0.65,Ireg_extrap[0]*0.9)
        chi_textitem.savex = q.max()*0.65
        chi_textitem.savey = Ireg_extrap[0]*0.9

        realRg_textitem = pg.TextItem("Rg_real : {0:7.2f}".format(sqrt(Rg_sq)), color='k', anchor=(0,0))
        realRg_textitem.setPos(r.max()*0.6,pr.max()*0.9)
        #realpars_textitem = pg.TextItem()

        self.prplot.clear()
        self.rawplot.clear()

        self.rawplot.addItem(self.Iq_errorbar)
        self.rawplot.plot(q,Iq,pen=self.redpen)
        self.rawplot.plot(q_full,Jreg_extrap,pen=self.blackpen)
        # only plot the unsmeared (Ireg) when data is smeared
        # otherwise the lines will just overlay and looks ugly
        if self.data.smeared:
            self.rawplot.plot(q_full,Ireg_extrap,pen=self.bluepen)
        self.rawplot.addItem(chi_textitem)

        self.prplot.addItem(realRg_textitem)
        self.prplot.plot(r,pr,pen=self.redpen)

        # after solving, plot the unsmeared normalized kratky plot
        # x = q*Rg
        # y = (q*Rg)^2 * I(q)/I(0)
        xc_kratky = q_full * self.data.Rg 
        yc_kratky = (xc_kratky)**2 * Ireg_extrap/Ireg_extrap[0]

        self.kratkyplot.clear()
        self.kratkyplot.plot(self.data.x_kratky, self.data.y_kratky, pen=self.redpen)
        self.kratkyplot.setLabel('bottom','q x Rg')
        self.kratkyplot.setLabel('left','(q*Rg)<sup>2</sup> x I(q)/I(0)')
        self.kratkyplot.addLine(x=1.73,pen=self.graydashedpen)
        self.kratkyplot.plot(xc_kratky, yc_kratky, pen=self.blackpen)        

        # save file (if user wants to)
        savefile = int(self.saveiftresult.checkState())

        if savefile>0:
            fn_out = str(self.iftresultout.text())
            
            # integrate P(r) to 1
            dr = r[1]-r[0]
            pr = pr/(pr.sum()*dr)

            fhd_out = open(fn_out,"wt")
            for n in range(Nr):
                fhd_out.write("{0:5.2f} {1:6.2E}\n".format(r[n],pr[n]))

            fhd_out.close()

            print "saved to {0}".format(fn_out)

    def biftgrideval(self):
        Nalpha = int(self.Nalpha.text())
        Ndmax  = int(self.Ndmax.text())
        minalpha = float(self.minlogalpha.value())
        maxalpha = float(self.maxlogalpha.value())
        mindmax  = float(self.mindmax.value())
        maxdmax  = float(self.maxdmax.value())
        alpharange = exp(linspace(minalpha,maxalpha,Nalpha))
        dmaxrange  = linspace(mindmax,maxdmax,Ndmax)
        Ncalc = len(alpharange) * len(dmaxrange)

        self.getiftparams()
        # run fortran wrapper ift from iftreg module
        Nr = self.data.Nr

        if self.data.smeared:
            y  = self.y
            Wy = self.Wy
        else:
            y = None
            Wy= None
            
        q  = self.data.q
        Iq = self.data.Iq
        sd = self.data.sd
        
        qmin_id = self.data.datarange[0]
        qmax_id = self.data.datarange[1]

        q  = q[qmin_id:qmax_id]
        Iq = Iq[qmin_id:qmax_id]
        sd = sd[qmin_id:qmax_id]

        self.statusbar.showMessage("Calculating Grid please wait ...")
        evimat = grideval(alpharange, dmaxrange, q, Iq, sd, Nr, y, Wy, self.weighdata, self.data.smeared)
        self.statusbar.showMessage("Grid done. Ready.")

        evidisp = exp(evimat - evimat.max())
        # get information about maxima in grid
        maxrow, maxcol = unravel_index(evidisp.argmax(), evidisp.shape)
        self.biftguess = [log(alpharange[maxrow]), dmaxrange[maxcol]]

        self.biftgridView.canvas.ax.cla()
        self.biftgridView.canvas.ax.contour(dmaxrange,log(alpharange),evidisp)
        self.biftgridView.canvas.ax.set_xlabel("Dmax")
        self.biftgridView.canvas.ax.set_ylabel("log(alpha)")
        self.biftgridView.canvas.draw()

    def biftsimplex(self):
        
        self.getiftparams()

        # run fortran wrapper ift from iftreg module
        Nr = self.data.Nr

        q  = self.data.q
        Iq = self.data.Iq
        sd = self.data.sd
        
        qmin_id = self.data.datarange[0]
        qmax_id = self.data.datarange[1]

        q  = q[qmin_id:qmax_id]
        Iq = Iq[qmin_id:qmax_id]
        sd = sd[qmin_id:qmax_id]

        smeared = self.smeared.checkState()

        p0 = self.biftguess

        if smeared>0:
            y = self.y
            Wy= self.Wy
        else:
            y = None
            Wy= None
            
        def objective(p0,q,Iq,sd,y,Wy,Nr,weighted):
            j1,j2,j3,j4,j5,j6,j7,evi = iftv2(exp(p0[0]),p0[1],q,Iq,sd,Nr,y,Wy,weighted,smeared)
            return -evi
        self.statusbar.showMessage("Optimizing Dmax & Alpha, please wait ...")
        res = minimize(objective, p0, args=(q,Iq,sd,y,Wy,Nr,self.weighdata), method="Nelder-Mead", tol=1e-2)
        self.statusbar.showMessage("Optimization done. Please check the 'Real Space' tab")

        optim_alpha = res.x[0]
        optim_Dmax  = res.x[1]

        self.logalpha_input.setText("{0:5.4f}".format(optim_alpha))
        self.dmax_input.setText("{0:5.2f}".format(optim_Dmax))


    def manualGuinier(self):

        qi = int(self.lb_guinier_spinbox.text())
        qf = int(self.ub_guinier_spinbox.text())
        self.data.manualGuinier(qi,qf)

        fit_x = self.data.guinierXY[0]
        fit_y = self.data.guinierXY[0]*self.data.guinierpars[0] + self.data.guinierpars[1]
        
        datlb = self.data.datarange[0]
        datub = self.data.datarange[1]

        guinier_x = self.data.q[datlb:datub]**2
        guinier_y = log(self.data.Iq[datlb:datub])

        self.data.kratky()

        self.guinierplot.clear()
        self.kratkyplot.clear()

        qRg_lower = self.data.guinierpts[0]
        qRg_upper = self.data.guinierpts[1]
        qRgfmt = "{0:5.2f} to {1:5.2f}".format(qRg_lower,qRg_upper)
        
        self.qRg_label.setText(qRgfmt)
        self.guinierRg_label.setText("{0:8.2f}".format(self.data.guinierRg))
        self.I0_label.setText("{0:8.4E}".format(self.data.guinierI0))

        self.guinierplot.plot(x=guinier_x, y=guinier_y,pen=self.redpen)
        self.guinierplot.setLabel('bottom','q<sup>2</sup>')
        self.guinierplot.setLabel('left','log[I(q)]')
        self.guinierplot.setXRange(0,fit_x.max()*1.4)
        self.guinierplot.setYRange(fit_y.min()*0.9, fit_y.max()*1.1)

        self.guinierplot.plot(x=fit_x,y=fit_y,pen=self.blackpen)
        self.kratkyplot.plot(self.data.x_kratky, self.data.y_kratky, pen=self.redpen)
        self.kratkyplot.setLabel('bottom','q x Rg')
        self.kratkyplot.setLabel('left','(q*Rg)<sup>2</sup> x I(q)/I(0)')
        self.kratkyplot.addLine(x=1.73,pen=self.graydashedpen)


def main():
    app = QApplication(sys.argv)
    app.setOrganizationName("Agard Lab @ UCSF")
    app.setApplicationName("SAXSgui v0.5")
    form = saxsgui_mainwindow()
    form.mainWindow.show()
    form.mainWindow.raise_()
    app.exec_()

# Execute main loop for GUI
if __name__ == "__main__":
    main()

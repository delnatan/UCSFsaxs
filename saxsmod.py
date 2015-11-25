from trans_smear import transc
from numpy import dot, vdot, array, log, where, exp, dot, outer, \
                pi, sin, linspace, eye, diag, ones, zeros, sqrt, abs, \
                vstack, hstack, repeat, arange, append
from numpy.linalg import inv, slogdet, norm
from sys import stdout
from scipy.optimize import minimize, nnls

# completely re-did the computation part, July 2nd, 2014. Daniel (Agard Lab)
# I figured that if you use newton update, iteration may converge in less than 20 steps!
# Solving the tikhonov problem is now implemented in python using numpy
# except when there's smearing correction, I use a fortran routine to speed up
# smearing matrix generation (for the grid-search/optimization)

# 07-20-2014 - the python routine in this module is used only for testing/debugging algorithm
#              the final algorithm is implemented in "IFTREG" module written in fortran


# internal use functions :
def fitline(x,y,w):
    # routine for fitting a straight line (x vs y) with with standard deviation
    # as weights, w.
    N = float(len(x))
    wi = 1.0/w**2 # inverse variance
    x_w = (wi*x).sum() / wi.sum()
    y_w = (wi*y).sum() / wi.sum()
    m = (wi*(x-x_w)*(y-y_w)).sum()/(wi*(x-x_w)**2).sum()
    b = y_w - m*x_w
    rmsd = sqrt( ((y-m*x-b)**2).sum() / N)
    return m,b,rmsd

class saxsdata:

    def __init__(self,filename,*args):
        self.filename = filename
        dat = open(filename,"rt")
        raw = dat.readlines()
        ext = filename[-3:].lower()
        
        if ext=="dat":
            try: #default format for reading data (used in beamlines, standard 3-column)
                self.q  = array([float(d.split()[0]) for d in raw if len(d.split())>1])
                self.Iq = array([float(d.split()[1]) for d in raw if len(d.split())>1])
                self.sd = array([float(d.split()[2]) for d in raw if len(d.split())>1])

            except:
                print("Possible header information, proceed to reading PDH...")
                Npts = int(raw[2].split()[0])
                self.q = array([float(d.split()[0]) for d in raw[5:(Npts+4)]])
                self.Iq= array([float(d.split()[1]) for d in raw[5:(Npts+4)]])
                self.sd= array([float(d.split()[2]) for d in raw[5:(Npts+4)]])

        elif ext=="pdh":
            Npts = int(raw[2].split()[0])
            self.q = array([float(d.split()[0]) for d in raw[5:(Npts+4)]])
            self.Iq= array([float(d.split()[1]) for d in raw[5:(Npts+4)]])
            self.sd= array([float(d.split()[2]) for d in raw[5:(Npts+4)]])

        elif ext=="the": # for theoretical calculations
            self.q  = array([float(d.split()[0]) for d in raw if len(d.split())>1])
            self.Iq = array([float(d.split()[1]) for d in raw if len(d.split())>1])
            self.sd = self.q * 0.0 + 1.0
        else:
            print "Data format {0} is not yet supported.".format(ext)

        self.check_and_convert_q_units()

        if len(args)==1:
            beamfn = args[0]
            beam = beamprofile(beamfn)

        # common properties of SAXS data    
        self.I0 = 0 # Real Space 
        self.Rg = 0 # Real Space
        self.r = 0
        self.pr = 0
        self.pr_error = 0
        self.guinierI0 = None
        self.guinierRg = None
        self.Ireg = None
        self.Jreg = None
        self.datarange = [0,len(self.q)]
        self.Nq_add = 0 # extrapolated scattering data length
        self.q_full = None

    def manualGuinier(self,qi,qf):
        qmin_id = self.datarange[0]
        qmax_id = self.datarange[1]

        q = self.q[qmin_id:qmax_id]
        Iq= self.Iq[qmin_id:qmax_id]
        sd= self.sd[qmin_id:qmax_id]

        qsq = q[qi:qf]**2
        LIq = log(Iq[qi:qf])
        Lw  = sd[qi:qf]/Iq[qi:qf]
        m,b,rmsd = fitline(qsq, LIq, Lw)
        Rg = sqrt(-3.0 * m)
        lower_qRg = q[qi]*Rg
        upper_qRg = q[qf]*Rg
        
        self.guinierRg = Rg
        self.guinierI0 = exp(b)
        self.guinierpts= [lower_qRg,upper_qRg]
        self.guinierpars = [m,b]
        self.guinierXY = [qsq,LIq]

    def kratky(self):
        # normalized kratky plot
        # x = q*Rg
        # y = (q*Rg)^2 * I(q)/I(0)
        # sqrt(3) ~ a peak for every well-folded protein
        qmin_id = self.datarange[0]
        qmax_id = self.datarange[1]
        y = (self.q[qmin_id:qmax_id] * self.guinierRg)**2 * self.Iq[qmin_id:qmax_id]/self.guinierI0
        x = (self.q[qmin_id:qmax_id] * self.guinierRg)
        self.y_kratky = y
        self.x_kratky = x

    def check_and_convert_q_units(self):

        qmax = self.q.max()
        if qmax<1.0: # then unit is in angstroms
            print "Unit is in angstroms.." 
        elif qmax>1.0:
            print "Unit is in nanometers .. multiplying by 0.1"
            self.q = 0.1*self.q

    def writeGNOM(self, savefilename):
        qi = self.datarange[0]
        qf = self.datarange[1]

        if self.solved:
            fout = open("%s_gnom.out" % savefilename,'wt')
            q,Iq,sd,r,pr,pr_error = self.q,self.Iq,self.sd,self.r,self.pr,self.pr_error
            if pr_error==0:
                pr_error = zeros(len(self.r))
            q  = q[qi:qf]
            Iq = Iq[qi:qf]
            sd = sd[qi:qf]
            q_full = self.q_full
            Jreg,Ireg = self.Jreg, self.Ireg
            Nq_add = len(q_full) - len(q)

            gheader = """
           ####    G N O M   ---   Version 4.6                       ####\n
   *******    Input file(s) : 

  Angular   range    :     from {0:9.4f} to {1:9.4f}
  Real space range   :     from {2:9.2f} to {3:9.2f}\n

      S          J EXP       ERROR       J REG       I REG\n\n"""

            fout.write(gheader.format(q.min(), q.max(),r.min(), r.max()))

            for i in range(0,Nq_add):
                fout.write("{0:12.4E}{1:36s}{2:12.4E}\n".format(q_full[i],' ',Ireg[i]))

            for i in range(Nq_add,len(Ireg)):
                fout.write("{0:12.4E}{1:12.4E}{2:12.4E}{3:12.4E}{4:12.4E}\n".format(\
                            q_full[i],Iq[i-Nq_add],sd[i-Nq_add],Jreg[i-Nq_add],Ireg[i]))

            prheader = """
           Distance distribution  function of particle  


       R          P(R)      ERROR\n\n"""

            fout.write(prheader)
            for i in range(len(r)):
                fout.write("{0:12.4E}{1:12.4E}{2:12.4E}\n".format(r[i],pr[i],pr_error[i]))


            footer = """          Reciprocal space: Rg = {0:8.2f}     , I(0) = {1:12.4E}
     Real space: Rg = {2:7.2f} +- {3:5.3f}  I(0) = {4:12.4E} +- {5:11.4E}\n"""

            fout.write(footer.format(self.guinierRg,self.guinierI0,self.Rg,0,self.Jreg[0],0))

            fout.close()
        else:   
            print "No solutions to be saved ... please solve the P(r) first.."

    def save_results(self, filename):
    	# things needed to save
    	# r,pr is final P(r) solution
    	# q,Iq,sd is experimental data (uncut), with data range (index used) for analysis in qi,qf
    	q,Iq,sd,r,pr = self.q,self.Iq,self.sd,self.r,self.pr
    	qi, qf = self.datarange
    	Jreg,Ireg = self.Jreg, self.Ireg

    	return False

class beamprofile:

    def __init__(self, filename):
        dat = open(filename, "rt")
        raw = dat.readlines()
        ext = filename[-3:].lower()

        if ext=="dat":
            Npts  = int(raw[2].split()[0])
            self.y = array([float(d.split()[0]) for d in raw[5:(Npts+4)]])
            self.Wy= array([float(d.split()[1]) for d in raw[5:(Npts+4)]])
            self.check_and_convert_q_units()
            self.dy = self.y[1] - self.y[0]
            self.normalize()
        elif ext=="txt":
            self.y = array([float(d.split()[0]) for d in raw])
            self.Wy= array([float(d.split()[1]) for d in raw])
            self.check_and_convert_q_units()
            self.dy = self.y[1] - self.y[0]
            self.normalize()
        else:
            print "Beam profile data: {0} is not recognized.".format(ext)

    def check_and_convert_q_units(self):
        ymax = self.y.max()
        if ymax>1.0:
            print "Beam profile unit is in nanometers .. multiplying by 0.1"
            self.y = 0.1*self.y

    def normalize(self):
        # normalize such that the integral is 0.5
        self.Wy = (self.Wy * 0.5) / (self.Wy.sum() * self.dy)

def trans(q, r, **kwargs):
    # make transformation matrix given for FT
    # the smearing construction is done following Steen Hansen's
    # Bayesian IFT program writteon in Fortran
    # I assume that the beam integral is equal to 0.5, and hence
    # the 2 constant at the last step of Matrix construction
    if len(kwargs) > 1:
        smear = True
        try:
            y  = kwargs['y']
            Wy = kwargs['Wy']
            dy = y[1]-y[0]
        except NameError:
            print "Please specify beam parameteres y and Wy"
            raise
    else:
        smear = False
    Nq = len(q)
    Nr = len(r)
    dr = r[1]-r[0]
    qr = outer(q,r)
    # sin(x)/x, when x=0, should be 1
    qr[qr<1.e-15] = 1.e-15
    A = 4*pi* sin(qr)/qr * dr   

    if smear:
    	# now this uses my fortran module 'transc' from 'trans_smear.f90'
        A = transc(q,r,y,Wy)
        # At = zeros((Nq,Nr))
        # for i in range(Nq):
        #     for j in range(Nr):
        #         qy = sqrt(q[i]**2 + y**2) * r[j]
        #         qy[qy<1e-15] = 1e-15
        #         wrk = Wy * sin(qy)/qy * dy
        #         At[i,j] = 2*dr* wrk.sum()
        # A = 4*pi*At
    return A

def spherepr(Nr, Dmax):
    # The integrated chord function of a sphere is 
    # Dmax**3/24 (V)
    # Since I(0) = 4*pi*V, we scale our P(r) so that the integral is 
    # equal to I(0)
    #V = (Iq0/(4*pi))*(24.0/Dmax**3)
    r = linspace(0,Dmax,Nr)
    dr = r[1]-r[0]
    chord = 1 - 1.5*(r/Dmax) + 0.5*(r/Dmax)**3
    pr = r**2 * chord
    pr = pr * (r<=Dmax)
    pr = pr/pr.sum()/dr
    return r,pr


def iftv2(alpha,Dmax,q,Iq,sd,Nr,y,Wy,weightdata,smeared):
    # calculates IFT solution and returns the complete set of solutions 
    # including extrapolation to zero angle
    Nq = len(Iq)
    
    # prepare q_full for extrapolation to zero angle
    dq = q[1]-q[0]
    q_extra = arange(0,q.min(),dq)
    q_full  = append(q_extra, q)

    if weightdata:
        sdnorm = 1./sd**2
        sdnorm = (sdnorm/sdnorm.sum()) * float(len(sdnorm))
        sdnorm = sdnorm.reshape(sdnorm.size,1)
    else:
        sdnorm = ones((Nq,1))

    r = linspace(0,Dmax,Nr)

    if smeared:
        K = transc(q,r,y,Wy)
        Kunsmeared = trans(q,r)
        Kfull = transc(q_full,r,y,Wy)
        Kfull_unsmeared = trans(q_full,r)
    else:
        K = trans(q,r)
        Kunsmeared = K
        Kfull = trans(q_full,r)
        Kfull_unsmeared = Kfull

    Kweighted = repeat(sdnorm,K.shape[1],axis=1) * K

    # construct matrix L and Z
    L = -0.5*eye(Nr,k=-1) + eye(Nr,k=0) - 0.5*eye(Nr,k=1)
    Z = zeros((Nr,Nr))
    Z[0,0] = 1.0
    Z[-1,-1] = 1.0

    C = vstack([Kweighted, alpha*L, 10.*alpha*Z])
    X = hstack([Iq*sdnorm[:,0], zeros(Nr), zeros(Nr)])

    sol,resnorm = nnls(C,X)
    pr = sol
    
    jreg = K.dot(sol)
    ireg = Kunsmeared.dot(sol)
    jreg_extrap = Kfull.dot(sol)
    ireg_extrap = Kfull_unsmeared.dot(sol)

    if weightdata:
        chisq = ((Iq-jreg)**2 / sd**2).mean()
    else:
        chisq = 1.0
        
    S0 = sum(-L.dot(sol)**2)
    U  = L + (K.T).dot(K)/alpha

    detsign,rlogdet = slogdet(U)
    rnorm = 0.5*(Nr*log(0.5)) + log(float(Nr+1))
    evidence = rnorm + (alpha*S0 - 0.5*chisq*Nq) - 0.5*rlogdet

    print "Chisq: {0:10.4f}\tEvidence:{1:10.5E}".format(chisq,evidence)

    return jreg, ireg, jreg_extrap, ireg_extrap, q_full, r, pr, evidence

def iftv3(K,Kunsmeared,alpha,Dmax,q,Iq,sd,Nr,weightdata,smeared):
    # third version of IFT in python, takes in two pre-calculated 
    # transformation matrices to speed up grid calculation
    # if unsmeared, do K=Ksmeared
    Nq = len(Iq)
    r = linspace(0,Dmax,Nr)
    # check if calculation is error-weighted
    # if no weight, multiply by diagonal of ones
    if weightdata:
        sdnorm = 1./sd**2
        sdnorm = (sdnorm/sdnorm.sum()) * float(len(sdnorm))
        sdnorm = sdnorm.reshape(sdnorm.size,1)
    else:
        sdnorm = ones((Nq,1))

    Kweighted = repeat(sdnorm,K.shape[1],axis=1) * K

    # construct matrix L and Z
    L = -0.5*eye(Nr,k=-1) + eye(Nr,k=0) - 0.5*eye(Nr,k=1)
    Z = zeros((Nr,Nr))
    Z[0,0] = 1.0
    Z[-1,-1] = 1.0

    C = vstack([Kweighted, alpha*L, 10.*alpha*Z])
    X = hstack([Iq*sdnorm[:,0], zeros(Nr), zeros(Nr)])

    sol,resnorm = nnls(C,X)
    pr = sol
    
    jreg = K.dot(sol)
    ireg = Kunsmeared.dot(sol)
    chisq = ((Iq-jreg)**2 / sd**2).mean()

    # Do calculation according to Hansen (Bayesian IFT)
    S0 = sum(-L.dot(sol)**2)
    U  = L + (K.T).dot(K)/alpha

    detsign,rlogdet = slogdet(U)
    rnorm = 0.5*(Nr*log(0.5)) + log(float(Nr+1))
    evidence = rnorm + (alpha*S0 - 0.5*chisq*Nq) - 0.5*rlogdet

    print "Chisq: {0:10.4f}\tEvidence:{1:10.5E}".format(chisq,evidence)

    return jreg, ireg, r, pr, evidence

def grideval(alpharange, dmaxrange, q, Iq, sd, Nr, y, Wy, weighdata, smeared):
    evigrid = zeros((len(alpharange),len(dmaxrange)))

    for j,d in enumerate(dmaxrange):
        r = linspace(0,d,Nr)
        if smeared:
            K = trans(q,r,y=y,Wy=Wy)
            Kunsmeared = trans(q,r)
        else:
            K = trans(q,r)
            Kunsmeared = K
        for i,a in enumerate(alpharange):
            #jreg,ireg,r,pr,evi = iftv2(a,d,q,Iq,sd,Nr,y,Wy,weighdata,smeared)
            jreg,ireg,r,pr,evi = iftv3(K,Kunsmeared,a,d,q,Iq,sd,Nr,weighdata,smeared)
            evigrid[i,j] = evi

    return evigrid

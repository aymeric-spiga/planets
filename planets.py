# modified by A. Spiga from files associated with R. Pierrehumbert's book
# added the possibility to use method associated with planets objects
# added the possibility to load a set of constants for a planet in a txt file
# added many functions and computations, and constants taken from scipy
# .... see tutorial

import os
import numpy as np
from datetime import datetime
import scipy.constants as cst

#Planetary database
#Source for planetary data, and some of the data on
#the moons, is http://nssdc.gsfc.nasa.gov/planetary/factsheet/

#-------------Basic physical constants-------------------
G = cst.G          # Gravitational constant
h = cst.h          # Planck's constant
sigma = cst.sigma  # Stefan-Boltzman constant
k = cst.k          # Boltzman thermodynamic constant
c = cst.c          # Speed of light
#-----------Thermodynamic constants----------------------
Avogadro = cst.Avogadro  # Avogadro's number
#Following will come out in J/(deg kmol), so
#that dividing Rstar by molecular weight gives
#gas constant appropriate for mks units
Rstarkilo = 1000.*k*Avogadro   #Universal gas constant
#-------------Useful planetary quantities----------------
astronomical_unit = cst.astronomical_unit       # astronomical unit in meters

desc = {}

#####################################
# convert from deg to rad
def deg_to_rad(angles): return angles*np.pi/180.

#####################################
# Planck function B_nu(T) or B_lambda(T)
# -- If lambda is given (in meters), the output units are W/(m^2 m)
def planck(temp,spec,kind="lambda"):
  if kind == "nu":
    X = (spec**3)/(c**2)
    Y = spec
  elif kind == "lambda":
    X = (c**2)/(spec**5)
    Y = c / spec
  print X
  print Y
  nnn = 2.*h*X
  ddd = np.exp(h*Y/(k*temp))-1.
  return nnn/ddd


#####################################
#####################################
#####################################
#####################################
class Planet:
    '''
    A Planet object contains basic planetary data.
   
    "print Planet.desc" for information

    For gas giants, "surface" quantities are given at the 1 bar level
    '''

############################################
### INIT
############################################

    #__repr__ object prints out a help string when help is
    #invoked on the planet object or the planet name is typed
    def __repr__(self):
        line1 =\
        'This planet object contains information on %s\n'%self.name
        line2 = 'Type \"help(Planet)\" for more information\n'
        return line1+line2

    def __init__(self):
        self.name = None ; desc["name"] = "Name of the planet"
        self.a = None ; desc["a"] = "Mean radius of planet (m)"
        self.mass = None ; desc["mass"] = "Mass of planet (kg)"
        self.g = None ; desc["g"] = "Surface gravitational acceleration (m/s**2)"
        self.L = None ; desc["L"] = "Annual mean solar constant (current) (W/m**2)"
        self.albedo = None ; desc["albedo"] = "Bond albedo (fraction)"
        self.rsm = None ; desc["rsm"] = "Semi-major axis of orbit about Sun (m)"
        self.year = None ; desc["year"] = "Sidereal length of year (s)"
        self.eccentricity = None ; desc["eccentricity"] = "Eccentricity (unitless)"
        self.day = None ; desc["day"] = "Mean tropical length of day (s)"
        self.obliquity = None ; desc["obliquity"] = "Obliquity to orbit (degrees)"
        self.Lequinox = None ; desc["Lequinox"] = "Longitude of equinox (degrees)"
        self.Tsbar = None ; desc["Tsbar"] = "Mean surface temperature (K)"
        self.Tsmax = None ; desc["Tsmax"] = "Maximum surface temperature (K)"
        self.M = None ; desc["M"] = "Molecular weight (g/mol)"
        self.cp = None ; desc["cp"] = "Specific heat capacity (J kg-1 K-1)"
        self.T0 = None ; desc["T0"] = "Typical atmospheric temperature (K)"
        self.incl = None ; desc["incl"] = "Orbit inclination (deg)"
        self.ascend = None ; desc["ascend"] = "Longitude of ascending node (deg)"
        self.omeg = None ; desc["omeg"] = "Argument of periapsis (deg)"
        self.date_peri = None ; desc["date_peri"] = "Date of perihelion"
        self.date_equi = None ; desc["date_equi"] = "Date of equinox"
        ## calculated
        self.R = None ; desc["R"] = "planetary gas constant"
        self.dryadiab = None ; desc["dryadiab"] = "dry adiabatic lapse rate"
        self.omega = None ; desc["omega"] = "planetary rotation rate"

############################################
### USEFUL METHODS FOR VALUES
############################################

    def show(self):
        # show objects attributes
        for k, v in vars(self).items():
          print k,v,desc[k]

    def convsecond(self):
        # convert earth days and hours in seconds
        if self.year is not None:
            self.year = self.year*24.*3600.
        if self.day is not None:
            self.day = self.day*3600.

    def convdate(self):
        # convert date peri and date equi in date format
        if self.date_peri is not None:
            self.date_peri = datetime.strptime(str(int(self.date_peri)), '%Y%m%d')
        if self.date_equi is not None:
            self.date_equi = datetime.strptime(str(int(self.date_equi)), '%Y%m%d')

    def calculate(self):
        # planetary gas constant
        if self.M is not None:
          self.R = Rstarkilo/self.M 
        else:
          self.R = None
        # adiabatic lapse rate
        if self.cp is not None:
          self.dryadiab = self.g/self.cp
        else:
          self.dryadiab = None
        # planetary rotation rate
        self.omega = 2.*np.pi/self.day

    def ini(self,name):
        # either have the file "name.txt" in /planet
        # ... or have it where you call
        string = "planets"
        whereset = "./"
        for path in os.environ['PYTHONPATH'].split(os.pathsep):
          if string in path: whereset = path + "/planet"
        if whereset[-1] != "/": whereset = whereset + "/"
        # be consistent
        self.name = name
        # set a dictionary with what's in the txt file
        cstplan = {}
        try: 
           f = open(whereset+name+".txt", 'r')
           for line in f:
            if "#" not in line and line != '\n' and line != '':
             variable,value = line.strip().split('=')
             cstplan[variable.strip()] = value.strip() 
           f.close()
        except IOError: 
           print "file not found: ",name+".txt" ; exit()
        # fill in object's attributes
        for k, v in vars(self).items():
          if k != "name":
            try:
              getval = cstplan[k]
              v = np.float(getval)
            except:
              #print k + " no value in file, set to None"
              v = None
            setattr(self,k,v)
        # do necessary converting
        self.convsecond()
        self.convdate()
        self.calculate()

############################################
### PHYSICAL CALCULATIONS as METHODS
############################################

    # Coriolis parameter
    def fcoriolis(self,lat=45.): 
        return 2.*self.omega*np.sin(deg_to_rad(lat))

    # calculate equivalent temperature
    def eqtemp(self,albedo=None):
        if albedo is None: albedo = self.albedo
        num = (1.-albedo)*self.L
        den = 4.*sigma
        return (num/den)**0.25

    # calculate Brunt-Vaisala frequency
    # ex: planets.Earth.N2(dTdz=-6.5e-3)
    # --> NB: dTdz could be an array
    def N2(self,T0=None,dTdz=None):
        if T0 is None: T0=self.T0
        if dTdz is None: dTdz=0.
        return (self.g / T0) * ( self.dryadiab + dTdz )

    # calculate scale height
    def H(self,T0=None):
        if T0 is None: T0=self.T0
        return self.R * T0 / self.g

    # calculate pseudo-altitude (log-pressure coordinates)
    def pseudoz(self,pressure,H=None,p0=None):
        if H is None: H=self.H()
        if p0 is None: p0=1.e5
        return H*np.log(p0/pressure)

    # planetary waves dispersion relationship
    def dispeqw(self,s,sigma,nu=0,lz=None,h=None,N2=None):
        a = self.a
        omega = self.omega
        g = self.g
        H = self.H()
        if N2 is None:
          N2 = self.N2()
        ##
        if h is None:
          if lz is None: 
            lz = H
          m = 2*np.pi/lz
          h = m**2 + (4*H*H)**(-1)
          h = g*h ; h = 1./h ; h = N2*h
        ##
        gamma = (4*a*a*omega*omega)/(g*h)
        lhs = (2.*nu+1.)*np.sqrt(gamma)
        w1 = np.where(sigma == 0.) ; sigma[w1] = np.nan
        func = gamma*(sigma**2) - s**2 - (s/sigma) - lhs
        return func
    
    # acosphi (pretty self-explanatory)
    # -- distance to axis of rotation
    def acosphi(self,lat):
        return self.a * np.cos(deg_to_rad(lat))

    # beta (Rossby parameter)
    # -- variation of f with latitude
    def beta(self,lat=None):
        if lat is None: lat=0.
        return 2 * self.omega * np.cos(deg_to_rad(lat)) / self.a

    # length of a degree of latitude / longitude (meters)
    # -- latitude if lat is not provided (or longitude@equator)
    # -- longitude if lat is provided
    def deglength(self,lat=None):
        if lat is None: lat=0.
        return np.cos(deg_to_rad(lat))*self.a*np.pi/180.

    # specific axial angular momentum
    # [= omega * a**2 if neither u, nor lat, is provided]
    def angmom(self,u=None,lat=None):
        if lat is None: lat=0.
        if u is None: u=0.
        acosphi = self.acosphi(lat)
        return acosphi*((self.omega*acosphi)+u)

    # local superotation index
    # -- excess local a.m. over u(equator) = 0
    # -- solid-body rotation: s = -1 (poles) to 0 (equator)
    # -- s>0 eddy forcing and transport of aam
    def superrot(self,u=None,lat=None):
        aam_norm = self.angmom(u=u,lat=lat) / self.angmom()
        return aam_norm-1.

    # Exner function
    # -- default p0 is 1 bar (as on Earth)
    # -- should be defining p0 in planetary constants?
    def exner(self,p,p0=1.e5):
        return (p/p0)**(self.R/self.cp)

    # potential temperature
    def tpot(self,temp,p,p0=1.e5):
        return temp/self.exner(p,p0=p0)

    # tanphi (pretty self-explanatory)
    def tanphi(self,lat):
        return np.tan(deg_to_rad(lat))

    # escape velocity
    # dePater & Lissauer 2.16
    def escape(self,r=None):
        if r is None: r = self.a
        return np.sqrt(2.*G*self.mass/r)

    # escape parameter 
    # dePater & Lissauer 4.76
    def lescape(self,r=None,m=None,temp=None):
        # defaut settings
        if r is None: 
          r = self.a # distance = planetary radius
        if temp is None: 
          temp = self.eqtemp(albedo=0.) # equ temp with A=0
        if m is None: 
          m = 1.007825e-3/Avogadro # atomic H
        # computations
        v0 = np.sqrt(2.*k*temp/m)
        ve = self.escape(r=r)
        return (ve/v0)**2

#----------------------------------------------------        
Earth = Planet() ; Earth.ini("Earth")
#----------------------------------------------------        
Mars = Planet() ; Mars.ini("Mars")
#----------------------------------------------------        
Saturn = Planet() ; Saturn.ini("Saturn")
#----------------------------------------------------
Venus = Planet() ; Venus.ini("Venus")
#----------------------------------------------------
Jupiter = Planet() ; Jupiter.ini("Jupiter")
#----------------------------------------------------
Pluto = Planet() ; Pluto.ini("Pluto")
#----------------------------------------------------
Mercury = Planet() ; Mercury.ini("Mercury")




#----------------------------------------------------        
Uranus = Planet()
Uranus.name = 'Uranus' #Name of the planet
Uranus.a = 25.362e6 #Mean radius of planet
Uranus.g = 8.87 #Surface gravitational acceleration
Uranus.albedo = .300 #Bond albedo
Uranus.L = 3.71 #Annual mean solar constant (current)
#
Uranus.rsm = 2872.46e9 #Semi-major axis
Uranus.year = 30685.4*24.*3600. #Sidereal length of year
Uranus.eccentricity = .0457 # Eccentricity
Uranus.day = 17.24*3600. #Mean tropical length of day
Uranus.obliquity = 97.77 #Obliquity to orbit (deg)
Uranus.Lequinox = None #Longitude of equinox (deg)
#
Uranus.Tsbar = 76. #Mean surface temperature
Uranus.Tsmax = None #Maximum surface temperature


#----------------------------------------------------        
Neptune = Planet()
Neptune.name = 'Neptune' #Name of the planet
Neptune.a = 26.624e6 #Mean radius of planet
Neptune.g = 11.15 #Surface gravitational acceleration
Neptune.albedo = .290 #Bond albedo
Neptune.L = 1.51 #Annual mean solar constant (current)
#
Neptune.rsm = 4495.06e9 #Semi-major axis
Neptune.year = 60189.0*24.*3600. #Sidereal length of year
Neptune.eccentricity = .0113 # Eccentricity
Neptune.day = 16.11*3600. #Mean tropical length of day
Neptune.obliquity = 28.32 #Obliquity to orbit (deg)
Neptune.Lequinox = None #Longitude of equinox (deg)
#
Neptune.Tsbar = 72. #Mean surface temperature
Neptune.Tsmax = None #Maximum surface temperature


#Selected moons

#----------------------------------------------------        
Moon = Planet()
Moon.name = 'Moon' #Name of the planet
Moon.a = 1.737e6 #Mean radius of planet
Moon.g = 1.62 #Surface gravitational acceleration
Moon.albedo = .11 #Bond albedo
Moon.L = 1367.6 #Annual mean solar constant (current)
#
Moon.rsm = Earth.rsm #Semi-major axis
Moon.year = Earth.year #Sidereal length of year
Moon.eccentricity = None # Eccentricity
Moon.day = 28.*24.*3600. #Mean tropical length of day (approx)
Moon.obliquity = None #Obliquity to orbit (deg)
Moon.Lequinox = None #Longitude of equinox (deg)
#
Moon.Tsbar = None #Mean surface temperature
Moon.Tsmax = 400. #Maximum surface temperature
Moon.Tsmin = 100. #Minimum surface temperature

Titan = Planet()
Titan.name = 'Titan' #Name of the planet
Titan.a = 2.575e6 #Mean radius of planet
Titan.g = 1.35 #Surface gravitational acceleration
Titan.L = Saturn.L #Annual mean solar constant (current)
Titan.albedo = .21 #Bond albedo (Not yet updated from Cassini)
#        
Titan.rsm = None #Semi-major axis
Titan.year = Saturn.year #Sidereal length of year
Titan.eccentricity = Saturn.eccentricity # Eccentricity ABOUT SUN
Titan.day = 15.9452*24.*3600. #Mean tropical length of day
Titan.obliquity = Saturn.obliquity #Obliquity to plane of Ecliptic
                                   #(Titan's rotation axis approx parallel
                                   # to Saturn's
Titan.Lequinox = Saturn.Lequinox #Longitude of equinox
#
Titan.Tsbar = 95. #Mean surface temperature
Titan.Tsmax = None #Maximum surface temperature

Europa = Planet()
Europa.name = 'Europa' #Name of the planet
Europa.a = 1.560e6 #Mean radius of planet
Europa.g = 1.31 #Surface gravitational acceleration
Europa.L = Jupiter.L #Annual mean solar constant (current)
Europa.albedo = .67 #Bond albedo
#        
Europa.rsm = Jupiter.rsm #Semi-major axis
Europa.year = Jupiter.year #Sidereal length of year
Europa.eccentricity = Jupiter.eccentricity # Eccentricity
Europa.day = 3.551*24.*3600. #Mean tropical length of day
Europa.obliquity = Jupiter.obliquity #Obliquity to plane of ecliptic
Europa.Lequinox = None #Longitude of equinox
#
Europa.Tsbar = 103. #Mean surface temperature
Europa.Tsmax = 125. #Maximum surface temperature

Triton = Planet()
Triton.name = 'Triton' #Name of the planet
Triton.a = 2.7068e6/2. #Mean radius of planet
Triton.g = .78 #Surface gravitational acceleration
Triton.L = Neptune.L #Annual mean solar constant (current)
Triton.albedo = .76 #Bond albedo
#        
Triton.rsm = Neptune.rsm #Semi-major axis
Triton.year = Neptune.year #Sidereal length of year
Triton.eccentricity = Neptune.eccentricity # Eccentricity about Sun
Triton.day = 5.877*24.*3600. #Mean tropical length of day
                             #Triton's rotation is retrograde
Triton.obliquity = 156. #Obliquity to ecliptic **ToDo: Check this.
                        #Note: Seasons are influenced by the inclination
                        #of Triton's orbit? (About 20 degrees to
                        #Neptune's equator
Triton.Lequinox = None #Longitude of equinox
#
Triton.Tsbar = 34.5 #Mean surface temperature
                    #This is probably a computed blackbody
                    #temperature, rather than an observation
Triton.Tsmax = None #Maximum surface temperature



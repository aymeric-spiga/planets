# modified by A. Spiga from files associated with R. Pierrehumbert's book
# added the possibility to use method associated with planets objects

import os
import numpy as np


#Planetary database
#Source for planetary data, and some of the data on
#the moons, is http://nssdc.gsfc.nasa.gov/planetary/factsheet/

#-------------Basic physical constants-------------------
#
#The following five are the most accurate 1986 values 
#
h = 6.626075540e-34    #Planck's constant
c = 2.99792458e8       #Speed of light
k = 1.38065812e-23      #Boltzman thermodynamic constant
sigma = 5.67051196e-8  #Stefan-Boltzman constant
G = 6.67428e-11        #Gravitational constant (2006 measurements)
#-----------Thermodynamic constants----------------------
#Following will come out in J/(deg kmol), so
#that dividing Rstar by molecular weight gives
#gas constant appropriate for mks units
N_avogadro = 6.022136736e23  #Avogadro's number
Rstarkilo = 1000.*k*N_avogadro   #Universal gas constant
#-------------Useful planetary quantities----------------
AU = 149597871000       # astronomical unit in meters

############################################################
desc = {}
desc["a"] = "Mean radius of planet (m)"
desc["g"] = "Surface gravitational acceleration (m/s**2)"
desc["L"] = "Annual mean solar constant (current) (W/m**2)"
desc["albedo"] = "Bond albedo (fraction)"
desc["rsm"] = "Semi-major axis of orbit about Sun (m)"
desc["year"] = "Sidereal length of year (s)"
desc["eccentricity"] = "Eccentricity (unitless)"
desc["day"] = "Mean tropical length of day (s)"
desc["obliquity"] = "Obliquity to orbit (degrees)"
desc["Lequinox"] = "Longitude of equinox (degrees)"
desc["Tsbar"] = "Mean surface temperature (K)"
desc["Tsmax"] = "Maximum surface temperature (K)"
desc["name"] = "Name of the planet"
desc["M"] = "Molecular weight (g/mol)"
desc["T0"] = "Typical atmospheric temperature (K)" 
desc["cp"] = "Specific heat capacity (J kg-1 K-1)"
desc["incl"] = "Orbit inclination (deg)"
desc["ascend"] = "Longitude of ascending node (deg)"
desc["omeg"] = "Argument of periapsis (deg)"
desc["date_peri"] = "Date of perihelion"
desc["date_equi"] = "Date of equinox"
#############################################################

class Planet:
    '''
    A Planet object contains basic planetary data.
   
    "print planets.desc" for information

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
        self.name = None #Name of the planet
        self.a = None #Mean radius of planet
        self.g = None #Surface gravitational acceleration
        self.L = None #Annual mean solar constant (current)
        self.albedo = None #Bond albedo
        
        self.rsm = None #Semi-major axis
        self.year = None #Sidereal length of year
        self.eccentricity = None # Eccentricity
        self.day = None #Mean tropical length of day
        self.obliquity = None #Obliquity to orbit
        self.Lequinox = None #Longitude of equinox

        self.Tsbar = None #Mean surface temperature
        self.Tsmax = None #Maximum surface temperature

        self.M = None #molar mass in grams per mol
        self.cp = None #specific heat capacity

        self.T0 = None #typical atm temperature

        self.incl = None #orbit inclination
        self.ascend = None #longitude ascending node
        self.omeg = None #argument of periapsis

        self.date_peri = None #date perihelion
        self.date_equi = None #date equinoxe

############################################
### USEFUL METHODS FOR VALUES
############################################

    def show(self):
        # show objects attributes
        for k, v in vars(self).items():
          print k,v,desc[k]

    def convsecond(self):
        # convert earth days and hours in seconds
        self.year = self.year*24.*3600.
        self.day = self.day*3600.

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

############################################
### PHYSICAL CALCULATIONS as METHODS
############################################

    def R(self): return Rstarkilo/self.M
        # planetary gas constant

    def dryadiab(self): return self.g/self.cp
        # adiabatic lapse rate

    def omega(self): return 2.*np.pi/self.day
        # planetary rotation rate

    def eqtemp(self): 
        # calculate equivalent temperature
        num = (1.-self.albedo)*self.L
        den = 4.*sigma
        return (num/den)**0.25

    def N2(self,T0=None,dTdz=None):
        # calculate Brunt-Vaisala frequency
        # ex: planets.Earth.N2(dTdz=-6.5e-3)
        # --> NB: dTdz could be an array
        if T0 is None: T0=self.T0
        if dTdz is None: dTdz=0.
        return (self.g / T0) * ( self.dryadiab() + dTdz )

    def H(self,T0=None):
        # calculate scale height
        if T0 is None: T0=self.T0
        out = self.R() * T0 / self.g
        return out

    def dispeqw(self,s,sigma,nu=0,lz=None,h=None,N2=None):
        a = self.a
        omega = self.omega()
        g = self.g
        H = self.H()
        if N2 is None:
          N2 = self.N2()
        ##
        if lz is None: 
          lz = H
        if h is None:
          m = 2*np.pi/lz
          h = m**2 + (4*H*H)**(-1)
          h = g*h ; h = 1./h ; h = N2*h
        ##
        gamma = (4*a*a*omega*omega)/(g*h)
        lhs = (2.*nu+1.)*np.sqrt(gamma)
        w1 = np.where(sigma == 0.) ; sigma[w1] = np.nan
        func = gamma*(sigma**2) - s**2 - (s/sigma) - lhs
        return func

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
Mercury = Planet()        
Mercury.name = 'Mercury' #Name of the planet
Mercury.a = 2.4397e6 #Mean radius of planet
Mercury.g = 3.70 #Surface gravitational acceleration
Mercury.albedo = .119 #Bond albedo
Mercury.L = 9126.6 #Annual mean solar constant (current)
#
Mercury.rsm = 57.91e9 #Semi-major axis
Mercury.year = 87.969*24.*3600. #Sidereal length of year
Mercury.eccentricity = .2056 # Eccentricity
Mercury.day = 4222.6*3600. #Mean tropical length of day
Mercury.obliquity = .01 #Obliquity to orbit (deg)
Mercury.Lequinox = None #Longitude of equinox (deg)
#
Mercury.Tsbar = 440. #Mean surface temperature
Mercury.Tsmax = 725. #Maximum surface temperature


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



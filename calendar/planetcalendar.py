#!/usr/bin/env python

from    numpy                 import    *
import  numpy                 as        np
import  matplotlib.pyplot     as        mpl
import  math
import  pylab
import datetime,time
from time import sleep
import sys
from datetime import datetime, timedelta
import planets

pi=math.pi

########################################################
def true_anomaly( M,e):
    
    # rad
    # initial approximation of eccentric anomaly
    E = M + e*sin(M)*(1.0 + e*cos(M))

    # iterate to improve accuracy
    E1 = E
    E = E1 - (E1 - e*sin(E1) - M)/(1 - e*cos(E1))
    eps=1.e-12
    while(abs(E-E1)>eps):
        E1 = E
        E = E1 - (E1 - e*sin(E1) - M)/(1 - e*cos(E1))
   
    # convert eccentric anomaly to true anomaly
    V = 2*arctan(pow((1 + e)/(1 - e),0.5)*tan(0.5*E))

    V = V%(2*pi) 
    
    return V


#######################################################
# INPUT :
#######################################################
def main(myplanet,date1,date2,day_step):
# date starting calendar
      #  date1=datetime(2000, 1, 1, 00, 0, 00)

# date ending calendar
      #  date2=datetime(2001, 2, 7, 0, 0, 0) 

# Calendar parameters
      #  day_step = 1   # step for calendar in Earth days

# PLOTS =1 'yes'
        plot=1

#######################################################
# ORBIT PARAMETERS:
#######################################################

	name=myplanet.name

	# date periapsis 
	date_peri=myplanet.date_peri  

	# date equinox
	date_equi=myplanet.date_equi  

	# orbital parameters
	e=myplanet.eccentricity   
	obl=myplanet.obliquity    
	#Mean tropical length of day (in seconds)
	day_sec=myplanet.day
	day_sol=day_sec/(24.*3600.)  # length of sol in earth days
	year_sec=myplanet.year # period in sec
	sma=myplanet.rsm/planets.AU

	JULIAN_EPOCH = datetime(2000, 1, 1, 12) # noon (the epoch is unrelated)
	J2000_JD = timedelta(2451545) # julian epoch in julian dates

	########################################################
	# Calcul True anomaly at equinox : vequi
	#######################################################

	# Calcul time between 1970 ref and periapsis (sec):
	timesecperi=time.mktime(date_peri.timetuple())
	# Calcul time between 1970 ref and time equinox (sec):
	timesecequi=time.mktime(date_equi.timetuple())

	# Calcul time between periapsis and equinox (sec):
	timeperiequi=timesecequi-timesecperi

	# calcul M at equinox since periapsis (deg):
	modequisecdiff=timeperiequi%year_sec
	Mequi=modequisecdiff/year_sec*360
	# M in radians
	Mequirad=Mequi*pi/180

	# Solving Kepler equation to get v at time 1 since periapsis (deg):
	vequi=true_anomaly(Mequirad,e)*180./pi

	########################################################
	# Calcul True anomaly at date 1 : v1
	#######################################################

	# Calcul time between 1970 ref and periapsis (sec):
	timesecperi=time.mktime(date_peri.timetuple())
	# Calcul time between 1970 ref and time 1 (sec):
	timesec1=time.mktime(date1.timetuple())
	# Calcul time between 1970 ref and time 2 (sec):
	timesec2=time.mktime(date2.timetuple())

	# Calcul time between periapsis and time 1 (sec):
	timesinceperi1=timesec1-timesecperi
	# Calcul time between time 1 and time 2 (sec):
	timesecdiff=timesec2-timesec1
	# Calcul number of step required:
	nbstep=int(timesecdiff/(day_step*24.*3600.))+1
	# calcul M at time 1 since periapsis (deg):
	modtimesecdiff=timesinceperi1%year_sec
	M1=modtimesecdiff/year_sec*360
	# M in radians
	Mrad=M1*pi/180

	# Solving Kepler equation to get v at time 1 since periapsis (deg):
	v1=true_anomaly(Mrad,e)*180./pi

	########################################################
	# TESTS INPUTS
	#######################################################

	if (timesec1>timesec2):
	   sys.exit("Wrong Input: Date 1 has to be lower than Date 2") 
	if e<0. or year_sec <0. or sma <0. or day_step<0. :
	   sys.exit("Wrong Input: negative parameters encountered") 
	if nbstep<1. :
	   sys.exit("Wrong Input: timestep too high compared to date 1 and 2 required") 

	########################################################

	print '  '
	print '***************************************'
	print '               Calendar                '
	print name.center(38, ' ')                    
	print '***************************************'
	print '  '
	print '--> date of periapsis : ',date_peri
	print '--> date of equinox : ',date_equi
	print '--> your choice of date for starting calendar : ',date1
	print '--> your choice of date for ending calendar : ',date2
	print '  '
	print ' ============================== '
	print '      Orbital parameters '
	print ' ============================== '
	print 'Eccentricity = ',e
	print 'Obliquity (deg) = ',obl
	print 'Year (in Earth year) = ',year_sec/(24.*3600.*365.25)
	print 'Period (Earth day) = ',year_sec/(24.*3600.)
	print 'Period (sec) = ',year_sec
	print 'True anomaly at equinox (deg) = ',vequi
	print '  '
	print ' ============================== '
	print '      Calendar parameters '
	print ' ============================== '
	print "Step of calendar (in Earth days) = ",day_step
	print "Number of steps = ",nbstep
	print "True anomaly at date 1 (deg) = ",v1
	print " "
	print " "

	#######################################################

	v=np.zeros(nbstep,dtype='f')
	ls=np.zeros(nbstep,dtype='f')
	M=np.zeros(nbstep,dtype='f')
	t=np.zeros(nbstep,dtype='f')
	tday=np.zeros(nbstep,dtype='f')
	tsinceperi=np.zeros(nbstep,dtype='f')
	decli=np.zeros(nbstep,dtype='f')
	had=np.zeros(nbstep,dtype='f')
	jj=np.zeros(nbstep,dtype='f')
	dist=np.zeros(nbstep,dtype='f')
	dates=np.zeros(nbstep,dtype='object')

	for i in range(nbstep):

	    tday[i]=i*day_step
	    t[i]=timesec1+tday[i]*24.*3600.
	    tsinceperi[i]=t[i]-timesecperi
	    # calcul M at time t since periapsis (deg):
	    M[i]=(tsinceperi[i]%year_sec)/year_sec*360.
	    # calcul v at time t since periapsis (deg):
	    v[i]=true_anomaly(M[i]*pi/180.,e)*180./pi  
	    ls[i]=(v[i]-vequi)%360.
	    dates[i]=datetime.fromtimestamp(t[i]).strftime('%Y-%m-%d %H:%M:%S')
	    decli[i]=arcsin(sin(obl*pi/180.)*sin(ls[i]*pi/180.))*180./pi
	    had[i]=tday[i]/day_sol
	    # Calcul Julian Date (Astronomical)
	    dt =datetime.strptime(dates[i],'%Y-%m-%d %H:%M:%S')
	    jj[i] = str((dt -  JULIAN_EPOCH + J2000_JD)).split()[0]
	    dist[i]=(sma*(1-e*e))/(1+e*cos(v[i]*pi/180.))
	#****** OUTPUT IN CALENDAR  *****

	fichier = open("calendar_"+name+".txt", "w")
	fichier.write(" time(day)     Sol            Earth date         JulianDate(As)    ls(deg)     M(deg)      v(deg)    decli(deg)     dist sol (UA)\n")
	fichier.write("********************************************************************************************************************************** \n")
	for j in range(nbstep):

        	fichier.write("%8s" % tday[j])
        	fichier.write('{0:>13}'.format(str('%9.3f'% had[j])))
        	fichier.write("%24s" % dates[j])
        	fichier.write('{0:>16}'.format(str('%14.3f'% jj[j])))
        	fichier.write('{0:>12}'.format(str('%7.3f'% ls[j])))
        	fichier.write('{0:>12}'.format(str('%7.3f'% M[j])))
        	fichier.write('{0:>12}'.format(str('%7.3f'% v[j])))
        	fichier.write('{0:>12}'.format(str('%7.3f'% decli[j])))
        	fichier.write('{0:>12}'.format(str('%7.3f'% dist[j])))
        	fichier.write("\n")
	fichier.write("\n")
	fichier.close()

        if (plot==1):
           mpl.figure()
           mpl.plot(had,decli)
           mpl.title('Declination Planet '+str(name))
           mpl.xlabel('Time (sol)',fontsize=20)
           mpl.ylabel('Declination (deg)',fontsize=20)
           mpl.grid(True)
           mpl.show()



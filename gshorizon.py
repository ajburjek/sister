import math

##Julien Date
#A FUNCTION WILL BE IMPLEMENTED
jd = 19255 #Oct 1 1995


##Local Sidereal Time (GMST)
ut = jd + 0.5
UT = ut - int(ut)
JD = jd - UT
TU = (JD - 2451545)/36525
GMST = 24110.54841+TU *(8640184.812866 + TU * (0.093104 - TU * 6.2E-6))
GMST = (GMST + 86400.0 * 1.00273790934 * UT) % 86400.0
ThetaG = 2*math.pi*GMST/86400.0 

print (ThetaG) #Radians

#Earth-Centered Inertial Reference Frame
#ADJUST LAT LONG
lat = 40* math.pi/180 #deg north
long = 75*math.pi/180 #deg west
alt = 0
Re = 6378.135
z = Re * math.sin(40) # z coordinate

R = (Re + alt) * math.cos(40)
theta = (ThetaG + long) % (2*math.pi)
x = R * math.cos(theta) #kilometers
y = R * math.sin(theta) #kilometers

print ('location of groundstation = ', x,y,z) #kilometers


#Topocentric-Horizon Coordinate System

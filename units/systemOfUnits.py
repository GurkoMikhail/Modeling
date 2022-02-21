# -*- C++ -*-
# ----------------------------------------------------------------------
# HEP coherent system of Units
#
# This file has been provided to CLHEP by Geant4 (simulation toolkit for HEP).
#
# The basic units are :
# millimeter              (millimeter)
# nanosecond              (nanosecond)
# Mega electron Volt      (MeV)
# positron charge         (eplus)
# degree Kelvin           (kelvin)
# the amount of substance (mole)
# luminous intensity      (candela)
# radian                  (radian)
# steradian               (steradian)
#
# Below is a non exhaustive list of derived and pratical units
# (i.e. mostly the SI units).
# You can add your own units.
#
# The SI numerical value of the positron charge is defined here,
# as it is needed for conversion factor : positron charge = e_SI (coulomb)
#
# The others physical constants are defined in the header file :
# PhysicalConstants.h
#
# Authors: M.Maire, S.Giani
#
# History:
#
# 06.02.96   Created.
# 28.03.96   Added miscellaneous constants.
# 05.12.97   E.Tcherniaev: Redefined pascal (to avoid warnings on WinNT)
# 20.05.98   names: meter, second, gram, radian, degree
#            (from Brian.Lasiuk@yale.edu (STAR)). Added luminous units.
# 05.08.98   angstrom, picobarn, microsecond, picosecond, petaelectronvolt
# 01.03.01   parsec    
# 31.01.06   kilogray, milligray, microgray    
# 29.04.08   use PDG 2006 value of e_SI
# 03.11.08   use PDG 2008 value of e_SI
# 19.08.15   added liter and its sub units (mma)
# 12.01.16   added symbols for microsecond (us) and picosecond (ps) (mma)
# 17.07.20   use PDG 2019 value of e_SI
# 06.05.21   added minute, hour, day, year and millielectronvolt (mma)

#
#
#
pi  = 3.14159265358979323846
twopi  = 2*pi
halfpi  = pi/2
pi2 = pi*pi

# 
# Length [L]
#
millimeter  = 1.
millimeter2 = millimeter*millimeter
millimeter3 = millimeter*millimeter*millimeter

centimeter  = 10.*millimeter   
centimeter2 = centimeter*centimeter
centimeter3 = centimeter*centimeter*centimeter
  
meter  = 1000.*millimeter  
meter2 = meter*meter
meter3 = meter*meter*meter

kilometer = 1000.*meter   
kilometer2 = kilometer*kilometer
kilometer3 = kilometer*kilometer*kilometer

parsec = 3.0856775807e+16*meter

micrometer = 1.e-6 *meter 
nanometer = 1.e-9 *meter
angstrom  = 1.e-10*meter
fermi = 1.e-15*meter

barn = 1.e-28*meter2
millibarn = 1.e-3 *barn
microbarn = 1.e-6 *barn
nanobarn = 1.e-9 *barn
picobarn = 1.e-12*barn

# symbols
nm  = nanometer
um  = micrometer

mm  = millimeter
mm2 = millimeter2
mm3 = millimeter3

cm  = centimeter   
cm2 = centimeter2
cm3 = centimeter3

liter = 1.e+3*cm3
L = liter
dL = 1.e-1*liter
cL = 1.e-2*liter
mL = 1.e-3*liter   

m  = meter  
m2 = meter2
m3 = meter3

km  = kilometer   
km2 = kilometer2
km3 = kilometer3

pc = parsec

#
# Angle
#
radian  = 1.  
milliradian = 1.e-3*radian
degree = (pi/180.0)*radian

steradian = 1.

# symbols
rad  = radian
mrad = milliradian
sr   = steradian
deg  = degree

#
# Time [T]
#
nanosecond  = 1.
second  = 1.e+9 *nanosecond
millisecond = 1.e-3 *second
microsecond = 1.e-6 *second
picosecond  = 1.e-12*second

minute = 60*second
hour   = 60*minute
day= 24*hour
year   = 365*day  

hertz = 1./second
kilohertz = 1.e+3*hertz
megahertz = 1.e+6*hertz

# symbols
ns = nanosecond
s = second
ms = millisecond
us = microsecond
ps = picosecond

#
# Electric charge [Q]
#
eplus = 1. # positron charge
e_SI  = 1.602176634e-19# positron charge in coulomb
coulomb = eplus/e_SI# coulomb = 6.24150 e+18 * eplus

#
# Energy [E]
#
megaelectronvolt = 1. 
electronvolt = 1.e-6*megaelectronvolt
kiloelectronvolt = 1.e-3*megaelectronvolt
gigaelectronvolt = 1.e+3*megaelectronvolt
teraelectronvolt = 1.e+6*megaelectronvolt
petaelectronvolt = 1.e+9*megaelectronvolt
millielectronvolt = 1.e-9*megaelectronvolt  

joule = electronvolt/e_SI# joule = 6.24150 e+12 * MeV

# symbols
MeV = megaelectronvolt
eV = electronvolt
keV = kiloelectronvolt
GeV = gigaelectronvolt
TeV = teraelectronvolt
PeV = petaelectronvolt

#
# Mass [E][T^2][L^-2]
#
kilogram = joule*second*second/(meter*meter)   
gram = 1.e-3*kilogram
milligram = 1.e-3*gram

# symbols
kg = kilogram
g = gram
mg = milligram

#
# Power [E][T^-1]
#
watt = joule/second# watt = 6.24150 e+3 * MeV/ns

#
# Force [E][L^-1]
#
newton = joule/meter# newton = 6.24150 e+9 * MeV/mm

#
# Pressure [E][L^-3]
#
#define pascal hep_pascal  # a trick to avoid warnings 
pascal = newton/m2   # pascal = 6.24150 e+3 * MeV/mm3
bar= 100000*pascal # bar= 6.24150 e+8 * MeV/mm3
atmosphere = 101325*pascal # atm= 6.32420 e+8 * MeV/mm3

#
# Electric current [Q][T^-1]
#
ampere = coulomb/second # ampere = 6.24150 e+9 * eplus/ns
milliampere = 1.e-3*ampere
microampere = 1.e-6*ampere
nanoampere = 1.e-9*ampere

#
# Electric potential [E][Q^-1]
#
megavolt = megaelectronvolt/eplus
kilovolt = 1.e-3*megavolt
volt = 1.e-6*megavolt

#
# Electric resistance [E][T][Q^-2]
#
ohm = volt/ampere# ohm = 1.60217e-16*(MeV/eplus)/(eplus/ns)

#
# Electric capacitance [Q^2][E^-1]
#
farad = coulomb/volt# farad = 6.24150e+24 * eplus/Megavolt
millifarad = 1.e-3*farad
microfarad = 1.e-6*farad
nanofarad = 1.e-9*farad
picofarad = 1.e-12*farad

#
# Magnetic Flux [T][E][Q^-1]
#
weber = volt*second# weber = 1000*megavolt*ns

#
# Magnetic Field [T][E][Q^-1][L^-2]
#
tesla = volt*second/meter2# tesla =0.001*megavolt*ns/mm2

gauss = 1.e-4*tesla
kilogauss = 1.e-1*tesla

#
# Inductance [T^2][E][Q^-2]
#
henry = weber/ampere# henry = 1.60217e-7*MeV*(ns/eplus)**2

#
# Temperature
#
kelvin = 1.

#
# Amount of substance
#
mole = 1.

#
# Activity [T^-1]
#
becquerel = 1./second 
curie = 3.7e+10 * becquerel
kilobecquerel = 1.e+3*becquerel
megabecquerel = 1.e+6*becquerel
gigabecquerel = 1.e+9*becquerel
millicurie = 1.e-3*curie
microcurie = 1.e-6*curie
Bq = becquerel
kBq = kilobecquerel
MBq = megabecquerel
GBq = gigabecquerel
Ci = curie
mCi = millicurie
uCi = microcurie

#
# Absorbed dose [L^2][T^-2]
#
gray = joule/kilogram 
kilogray = 1.e+3*gray
milligray = 1.e-3*gray
microgray = 1.e-6*gray

#
# Luminous intensity [I]
#
candela = 1.

#
# Luminous flux [I]
#
lumen = candela*steradian

#
# Illuminance [I][L^-2]
#
lux = lumen/meter2

#
# Miscellaneous
#
perCent = 0.01 
perThousand = 0.001
perMillion  = 0.000001


# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 13:38:28 2021

@author: Louis Huchet

"""

# # IP 2020-21: DEVELOPING A VPP FOR A NACRA 17

# ## 1) Aerodynamic module

# Import necssary packages to run the code:

# In[2]:


import math
import numpy as np
import pandas as pd
import time
from math import degrees as deg
from math import radians as rad
from math import cos as cos
from math import sin as sin
from math import atan as atan
from math import sqrt as sqrt
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import scipy.interpolate as spi
import matplotlib.pyplot as plt
from matplotlib.pyplot import polar
from IPython.display import display
pd.set_option('display.max_rows', 40)
pd.set_option('display.max_columns', 20)


# ### Initial data (boat, environment...)

# In[3]:


LWL = 5.25 #m
rho_water = 1025 #Kg/m3
rho_air = 1.225 #Kg/m3
A_main = 14.45 #m²
A_jib = 4 #m²
A_upwind = A_main + A_jib
A_spi = 18.5 #m²
A_downwind = A_main + A_jib + A_spi
AR = 4.85 #aspect ratio of mainsail
k = 1/(math.pi*AR)
nu_air = 1.802e-5 #kg/m-s
nu_water = 1.1892e-6
v = 1.48e-5 #m²/s
RM_max = 7397.24 #N.m
PM_max = 4550 #N.m
g = 9.81 #kg/s²
boat_weight = 163*g
crew_weight = 120*g
hull_form = 1.22
Aw_1hull = 1.914 #m²
Aw_2hulls = 3.828 #m²


# ### Aerodynamic sail coefficients, based on ORC VPP data

# In[4]:


Cl_main = [0,0.86207,1.05172,1.16379,1.34698,1.35345,1.26724,0.93103,0.38793,-0.11207]
Cd_main = [0.0431,0.02586,0.02328,0.02328,0.03259,0.11302,0.3825,0.96888,1.31578,1.34483]
Beta_main = [0,7,9,12,28,60,90,120,150,180]
a = np.array([Beta_main, Cl_main, Cd_main])
Cl_jib = [0,1,1.375,1.45,1.45,1.25,0.4,0,-0.1]
Cd_jib = [0.05,0.032,0.031,0.037,0.25,0.35,0.73,0.95,0.9]
Beta_jib = [7,15,20,27,50,60,100,150,180]
b = np.array([Beta_jib, Cl_jib, Cd_jib])
df_spi = pd.read_excel("Copy of Database.xlsx", "Asymmetric spinnaker on pole", engine = 'openpyxl')
Cl_spi = np.array(df_spi['Cl'])
Cd_spi = np.array(df_spi['Cd'])
Beta_spi = np.array(df_spi['beta'])


# ### Interpolated values for a greater set ofAWA (Main, jib, Spinnaker)

# In[5]:


beta_vals = np.linspace(0,180,40) ##increment of 4.5°
z = list(beta_vals) +Beta_main + Beta_jib
set(z)
z.sort()
z = list(dict.fromkeys(z))
beta_vals_spi = list(dict.fromkeys(list(np.linspace(28, 180, 40))+list(Beta_spi)))
beta_vals_spi.sort()

Cl_main_interp = np.interp(z, Beta_main, Cl_main)
Cd_main_interp = np.interp(z, Beta_main, Cd_main)
Cl_jib_interp = np.interp(z, Beta_jib, Cl_jib)
Cd_jib_interp = np.interp(z, Beta_jib, Cd_jib)
Cl_spi_interp = np.interp(beta_vals_spi, Beta_spi, Cl_spi)
Cd_spi_interp = np.interp(beta_vals_spi, Beta_spi, Cd_spi)

#plt.plot(z, Cl_main_interp, '-k', marker = 'x', markevery = 5, markersize = 4, label = '$C_L$ Mainsail')
#plt.plot(z, Cd_main_interp, '--k', marker = 'x', markevery = 5, markersize = 4, label = '$C_D$ Mainsail')
#plt.plot(z, Cl_jib_interp, '-k', marker = '^', markevery = 5, markersize = 4, label = '$C_L$ Jib')
#plt.plot(z, Cd_jib_interp, '--k', marker = '^', markevery = 5, markersize = 4, label= '$C_D$ Jib')
#plt.plot(beta_vals_spi, Cl_spi_interp, '-k', marker = 's', markevery = 5, markersize = 4, label = '$C_L$ Spinnaker')
#plt.plot(beta_vals_spi, Cd_spi_interp, '--k', marker = 's', markevery = 5, markersize = 4, label = '$C_D$ Spinnaker')
#plt.xlabel('Apparent wind angle')
#plt.ylabel('Coefficient')
#plt.legend(bbox_to_anchor=(1, 1))

#
#plt.subplot(3, 2, 1)
#plt.plot(z, Cl_main_interp, 'x')
#plt.title("Cl mainsail vs. wind angle")
#plt.subplot(3, 2, 2)
#plt.plot(z, Cd_main_interp, 'x')
#plt.title("Cd mainsail vs. wind angle")
#plt.subplot(3, 2, 3)
#plt.plot(z, Cl_jib_interp, 'x')
#plt.title("Cl jib vs. wind angle")
#plt.subplot(3, 2, 4)
#plt.plot(z, Cd_jib_interp, 'x')
#plt.title("Cd jib vs. wind angle")
#plt.subplot(3,2,5)
#plt.plot(beta_vals_spi, Cl_spi_interp, 'x')
#plt.title("Cl spinnaker vs wind angle")
#plt.subplot(3,2,6)
#plt.plot(beta_vals_spi, Cd_spi_interp, 'x')
#plt.title("Cd spinnaker vs. wind angle")
#plt.rcParams["figure.figsize"] = (15,10)
#plt.tight_layout()
#plt.show()


# ### Combined Interpolated values for main+jib sails

# In[6]:


Cl_upwind = [x*A_main/A_upwind + y*A_jib/A_upwind for x, y in zip(Cl_main_interp, Cl_jib_interp)]
Cl_upwind[0] = 0.0001 #avoids math error
Cdp_upwind = [x*A_main/A_upwind + y*A_jib/A_upwind for x, y in zip(Cd_main_interp, Cd_jib_interp)]
Cdi = [x**2/(AR*math.pi) for x in Cl_upwind]
Cd_upwind = [x + y for x, y in zip(Cdp_upwind, Cdi)]
#
#plt.figure(1)
#plt.plot(z, Cl_upwind, 'x', label="Cl")
#plt.title("Combined Main + jib Cl & Cd coeffs vs. AWA")
#plt.plot(z, Cd_upwind, '^', label="Cd")
#plt.rcParams["figure.figsize"] = (12,5)
#plt.tight_layout()
#plt.show()
#
#
##In[7]:
#
#plt.figure(4)
#plt.plot(Cd_upwind, Cl_upwind, '-x')
#plt.title("Cl vs Cd (Main+jib combined)")
#plt.xlabel("Cd")
#plt.ylabel("Cl")
#plt.rcParams["figure.figsize"] = (6,5)
#plt.show()


# ### Combined interpolated values for Main + jib + Spinnaker (induced drag not included)

# In[8]:


Cl_downwind = [(x*A_main+y*A_jib+z*A_spi)/A_downwind for x, y, z in zip(Cl_main_interp, Cl_jib_interp, Cl_spi_interp)]
Cd_downwind = [(x*A_main+ y*A_jib+z*A_spi)/A_downwind for x, y, z in zip(Cd_main_interp, Cd_jib_interp, Cd_spi_interp)]
Cdi_downwind = [x**2/(AR*math.pi) for x in Cl_downwind]
Cd_downwind = [x + y for x, y in zip(Cd_downwind, Cdi_downwind)]

#plt.figure(1)
#plt.plot(beta_vals_spi, Cl_downwind, 'x', label="Cl Spi")
#plt.title("Combined Main + jib + spinnaker Cl coeff. vs wind angles")
#plt.plot(beta_vals_spi, Cd_downwind, '^', label="Cd Spi")

#plt.figure(3)
#plt.plot(Cd_downwind, Cl_downwind, 'x')
#plt.title("Cl vs Cd (Main+jib+spinnaker combined)")
#plt.xlabel("Cd")
#plt.ylabel("Cl")
#plt.rcParams["figure.figsize"] = (12,5)
#plt.tight_layout()
#plt.show()


# ## 2) Maxsurf temporary data

# ### Hull resistance dataframe and Stability 

# In[12]:


df = pd.read_excel("Copy of Database.xlsx", "120kg crew hull resistance", engine = 'openpyxl')
speed = np.array(df['Speed (m/s)'])


def get_hull_R(Vs, crew_weight, heel_angle):
    if crew_weight == 120*g and heel_angle > 4:
        return np.interp(Vs, speed, np.array(pd.read_excel("Copy of Database.xlsx", "120kg crew hull resistance", engine = 'openpyxl')['Slender body resistance 1 hull (N)'], dtype = np.float))
    elif crew_weight == 120*g and heel_angle <= 4:
        return np.interp(Vs, speed, np.array(pd.read_excel("Copy of Database.xlsx", "120kg crew hull resistance", engine = 'openpyxl')['Slender body resistance 2 hulls (N)'], dtype = np.float))
    elif crew_weight == 150*g and heel_angle > 4:
        return np.interp(Vs, speed, np.array(pd.read_excel("Copy of Database.xlsx", "150kg crew hull resistance", engine = 'openpyxl')['Slender body resistance 1 hull (N)'], dtype = np.float))
    elif crew_weight == 150*g and heel_angle <= 4:
        return np.interp(Vs, speed, np.array(pd.read_excel("Copy of Database.xlsx", "150kg crew hull resistance", engine = 'openpyxl')['Slender body resistance 2 hulls (N)'], dtype = np.float))
    elif crew_weight == 180*g and heel_angle > 4:
        return np.interp(Vs, speed, np.array(pd.read_excel("Copy of Database.xlsx", "180kg crew hull resistance", engine = 'openpyxl')['Slender body resistance 1 hull (N)'], dtype = np.float))
    elif crew_weight == 180*g and heel_angle <= 4:
        return np.interp(Vs, speed, np.array(pd.read_excel("Copy of Database.xlsx", "180kg crew hull resistance", engine = 'openpyxl')['Slender body resistance 2 hulls (N)'], dtype = np.float))


df_RM = pd.read_excel("Copy of Database.xlsx", "Stability case 120kg crew", engine = 'openpyxl')
heel = np.array(df_RM['Heel (deg)'])[:15]
Fn_list = np.array(pd.read_excel("Copy of Database.xlsx", "Resid Resistance Molland-Soton", engine = 'openpyxl')['Fn'], dtype = np.float) 
Rr_1hull = np.array(pd.read_excel("Copy of Database.xlsx", "Resid Resistance Molland-Soton", engine = 'openpyxl')['Total Rr 1hull [N]'], dtype = np.float)
Rr_2hulls = np.array(pd.read_excel("Copy of Database.xlsx", "Resid Resistance Molland-Soton", engine = 'openpyxl')['Total Rr 2hulls [N]'], dtype=np.float)


# In[15]:


#RM_tot_0trap = np.array(df_RM['Total RM 0 trap (N.m)'])
#RM_tot_1trap = np.array(df_RM['Total RM 1 trap (N.m)'])
#RM_tot_2trap = np.array(df_RM['Total RM 2 trap (N.m)'])
heel_gz = np.array(pd.read_excel("Copy of Database.xlsx", "Stability case 120kg crew", engine = 'openpyxl')['Heel (deg)'])[:15]
GZ_list = np.array(pd.read_excel("Copy of Database.xlsx", "Stability case 120kg crew", engine = 'openpyxl')['GZ (m)'])

if crew_weight == 120*g: 
    RM_tot_0trap = np.array(pd.read_excel("Copy of Database.xlsx", "Stability case 120kg crew", engine = 'openpyxl')['Total RM 0 trap (N.m)'], dtype = np.float)[:15]  #N.m
    RM_tot_1trap = np.array(pd.read_excel("Copy of Database.xlsx", "Stability case 120kg crew", engine = 'openpyxl')['Total RM 1 trap (N.m)'], dtype = np.float)[:15]  #N
    RM_tot_2trap = np.array(pd.read_excel("Copy of Database.xlsx", "Stability case 120kg crew", engine = 'openpyxl')['Total RM 2 trap (N.m)'], dtype = np.float)[:15]  #N
    RM_max = max(np.array(pd.read_excel("Copy of Database.xlsx", "Stability case 120kg crew", engine = 'openpyxl')['Total RM 2 trap (N.m)'], dtype = np.float)) #N.m
if crew_weight == 150*g:
    RM_tot_0trap = np.array(pd.read_excel("Copy of Database.xlsx", "Stability case 150kg crew", engine = 'openpyxl')['Total RM 0 trap (N.m)'], dtype = np.float)[:15]  #N.m
    RM_tot_1trap = np.array(pd.read_excel("Copy of Database.xlsx", "Stability case 150kg crew", engine = 'openpyxl')['Total RM 1 trap (N.m)'], dtype = np.float)[:15]  #N
    RM_tot_2trap = np.array(pd.read_excel("Copy of Database.xlsx", "Stability case 150kg crew", engine = 'openpyxl')['Total RM 2 trap (N.m)'], dtype = np.float)[:15]  #N
    RM_max = max(np.array(pd.read_excel("Copy of Database.xlsx", "Stability case 150kg crew", engine = 'openpyxl')['Total RM 2 trap (N.m)'], dtype = np.float))  #N.m
if crew_weight == 180*g:
    RM_tot_0trap = np.array(pd.read_excel("Copy of Database.xlsx", "Stability case 180kg crew", engine = 'openpyxl')['Total RM 0 trap (N.m)'], dtype = np.float)[:15]  #N.m
    RM_tot_1trap = np.array(pd.read_excel("Copy of Database.xlsx", "Stability case 180kg crew", engine = 'openpyxl')['Total RM 1 trap (N.m)'], dtype = np.float)[:15]  #N
    RM_tot_2trap = np.array(pd.read_excel("Copy of Database.xlsx", "Stability case 180kg crew", engine = 'openpyxl')['Total RM 2 trap (N.m)'], dtype = np.float)[:15]  #N
    RM_max = max(np.array(pd.read_excel("Copy of Database.xlsx", "Stability case 180kg crew",engine = 'openpyxl')['Total RM 2 trap (N.m)'], dtype = np.float))  #N.m





# In[16]:


#ax1 = plt.subplot()
#ax1.set_xlabel('heel angle (deg)')
#ax1.set_ylabel('Righting moment [N.m]')
#ax1.plot(heel_gz, RM_tot_0trap, '-r')
#ax1.plot(heel_gz, RM_tot_1trap, '-g')
#ax1.plot(heel_gz, RM_tot_2trap, '-y')
#ax1.legend(['0 trap', '1 trap', '2 trap'])
#
#ax2 = ax1.twinx()
#ax2.set_ylabel('GZ [m]')
#ax2.plot(heel_gz, GZ_list, 'bx')
#ax2.tick_params(axis='y', labelcolor='blue')
#plt.title('Righting moment vs heel angle')
#plt.tight_layout()
#plt.show()


# ### List of intermediate functions

# 1) Determine the height of the aerodynamic centre of effort. Function of previous CE height, flattening factor and heel angle. 

# In[17]:


def getZce(Zce, flat):
    return Zce*(1-0.203*(1-flat) - 0.451*(1-flat))


# In[19]:


def getAWA(Vt, gamma, Vs):
    awa = deg(atan(sin(rad(gamma))/((Vs/Vt)+cos(rad(gamma)))))
    if awa < 0:
        return 180+awa
    else:        
        return awa   #result in degrees

def getAWS(Vt, gamma, Vs):
    beta = getAWA(Vt, gamma, Vs)
    V_app = Vt*sin(rad(gamma))/sin(rad(beta))
    return V_app


# 5) Determine the aerodynamic sail lift and drag coefficients of combined sails for a given AWA.

# In[20]:


def getCl(beta_deg):
    return np.interp(beta_deg, z, Cl_upwind)

def getCd(beta_deg):
    return np.interp(beta_deg, z, Cd_upwind)

def getCl_downwind(beta_deg):
    return np.interp(beta_deg, beta_vals_spi, Cl_downwind)

def getCd_downwind(beta_deg):
    return np.interp(beta_deg, beta_vals_spi, Cd_downwind)

# 6) Determine aerodynamic Lift and Drag forces generated from the sails. Function of Apparent wind speed, angle and sail area. Assumes optimised sail trim 

# In[21]:


def getLift(Va, beta_deg, gamma):
    if gamma < 90:
#        if (Va < 10 and beta_deg < 40) or (Va >= 10 and beta_deg < 75):
        Cl = getCl(beta_deg)
        return 0.5*rho_air*Cl*A_upwind*(Va)**2
#        elif (Va < 10 and beta_deg > 40) or (Va >= 10 and beta_deg > 75):
#        Cl = getCl_downwind(beta_deg)
#            return 0.5*rho_air*Cl*A_downwind*(Va)**2
    else:
        Cl = getCl_downwind(beta_deg)
        return 0.5*rho_air*Cl*A_downwind*(Va)**2

    
def getDrag(Va, beta_deg, gamma):
    if gamma < 90:
#        if (Va < 10 and beta_deg < 40) or (Va >= 10 and beta_deg < 75):
        Cd = getCd(beta_deg)
        return 0.5*rho_air*Cd*A_upwind*(Va)**2
#        elif (Va < 10 and beta_deg > 40) or (Va >= 10 and beta_deg >75):
#            Cd = getCd_downwind(beta_deg)
#        return 0.5*rho_air*Cd*A_downwind*(Va)**2
    else:
        Cd = getCd_downwind(beta_deg)
        return 0.5*rho_air*Cd*A_downwind*(Va)**2


# 7) Determine aerodynamic driving and heeling forces, Fx and Fy, generated by the sails. Function of the lift and drag generated by the sails

# In[22]:


def getFx(Va, beta_deg, gamma, flat_fact):
    L_sails = flat_fact*getLift(Va, beta_deg, gamma)
    D_sails = getDrag(Va, beta_deg, gamma)
    beta_rad = rad(beta_deg)
    return L_sails*sin(beta_rad) - D_sails*cos(beta_rad)

def getFy(Va, beta_deg, gamma, flat_fact):
    L_sails = flat_fact*getLift(Va, beta_deg, gamma)
    D_sails = getDrag(Va, beta_deg, gamma)
    beta_rad = math.radians(beta_deg)
    return L_sails*cos(beta_rad) + D_sails*sin(beta_rad)


# 8) Determine optimised angle of attack of the sails for a given AWA (from aerodynamic coefficients polar Cl/Cd).

# In[23]:


def getAoA(beta_deg):
    return deg(atan(getCd(beta_deg)/getCl(beta_deg)))


# 9) Determine Heeling moment. function of heeling force and aerodynamic CE height:

# In[24]:


def getHM1(Fy, Zce):
    return Fy*(Zce+0.3)


# 10) Determine heel angle from righting moment dataframe (based on GZ curve and righting moment created by number of crew out on trapeze (0, 1 or 2). Assume equilibrium condition RM = HM if HM is smaller than the maximum righting moment)

# In[25]:


def get_heelangle(HM_1):
    if crew_weight == 120*g:
        if 0 < HM_1 < 4200: #0 crew on trapeze
            return np.interp(HM_1, RM_tot_0trap, heel)
        elif 4200 <= HM_1 <= 4950: #1 crew on trapeze
            return np.interp(HM_1, RM_tot_1trap, heel)
        elif 4950 < HM_1 <= RM_max: #2 crew on trapeze
            return np.interp(HM_1, RM_tot_2trap, heel)
    if crew_weight == 150*g:
        if 0 < HM_1 < 4900: #0 crew on trapeze
            return np.interp(HM_1, RM_tot_0trap, heel)
        elif 4900 <= HM_1 <= 5750: #1 crew on trapeze
            return np.interp(HM_1, RM_tot_1trap, heel)
        elif 5750 < HM_1 <= RM_max: #2 crew on trapeze
            return np.interp(HM_1, RM_tot_2trap, heel)
    if crew_weight == 180*g:
        if 0 < HM_1 < 5530: #0 crew on trapeze
            return np.interp(HM_1, RM_tot_0trap, heel)
        elif 5530 <= HM_1 <= 6575: #1 crew on trapeze
            return np.interp(HM_1, RM_tot_1trap, heel)
        elif 6575 < HM_1 <= RM_max: #2 crew on trapeze
            return np.interp(HM_1, RM_tot_2trap, heel)       

# 11) Determine frictionnal hull resistance

# In[58]:


def getRf(Vs, heel_angle):
    Re = Vs*0.7*LWL/nu_water
    Cf = 0.075/(math.log10(Re)-2)**2
    if heel_angle < 4:
        return hull_form*(Cf*0.5*rho_water*Aw_2hulls*(Vs)**2)
    else:
        return hull_form*(Cf*0.5*rho_water*Aw_1hull*(Vs)**2)

# 12) Determine Residuary Resistance

# In[59]:


def getRr(Vs, heel_angle):
    Fn = Vs/sqrt(g*LWL)
    if heel_angle < 4:
        return np.interp(Fn, Fn_list, Rr_2hulls)
    else:
        return np.interp(Fn, Fn_list, Rr_1hull)

# Determine Heel resistance

# In[60]:

def getheelR(Rh, heel_angle):
    if heel_angle >4:
        return 0.05*Rh
    else: 
        return 0

# 13) Determine Total Canoe body Resistance 

# In[61]:


def getRh(Rf, Rr):
    return Rf+Rr

# 13) Determine the appendage resistance + windage. Assume now only 7% of Total hull resistance for ease and 5% for windage

# In[62]:


#def get_appendageR(Rh):
#    return 0.1*Rh


# In[63]:


def getwindage_hull(Va, beta):
    W=0.39 #Dimensions of the hull Width, height and length
    H=0.375
    l = 5.25
    hull_area = (2*cos(rad(beta))*W + l*sin(rad(beta)))*H
    Cdx = 0.4
    Cdy = 0.9
    Cd_hull = Cdx*cos(rad(beta))+Cdy*sin(rad(beta))
    return cos(rad(beta))*0.5*rho_air*(Va)**2*Cd_hull*hull_area

def get_total_windage(Va, beta):
    windage_hull = getwindage_hull(Va, beta)
    Cd_tramp = 0.05
    Cd_crew = 1.02
    Cd_rigging = 1.2
    A_crew = 1.6
    A_rigging = 0.116
    A_tramp = 2.7
    windage_tramp = cos(rad(beta))*0.5*rho_air*(Va)**2*Cd_tramp*A_tramp
    windage_crew = cos(rad(beta))*0.5*rho_air*(Va)**2*Cd_crew*A_crew
    windage_rigging = cos(rad(beta))*0.5*rho_air*(Va)**2*Cd_rigging*A_rigging
    return windage_hull + windage_tramp + windage_crew + windage_rigging


def getPM(Fx, Zce):
    return Fx*Zce


# In[65]:

"""APPENDAGE DRAG"""
c_dagg = 0.238 #chord daggerboard's foils
c_rudd = 0.200 #chord rudder
t_dagg = 0.03808 #max thickness daggerboards
t_rudder = 0.12*c_rudd #assume NACA0012 profile (thickness = 12% of chord)
alpha_0_rudder = -1.8
s_rudder = 0.5
A_dagg_up = 0.284 #Total area of one surface facing upward on 1 daggerboard from baseline to tip of the foil
A_tot_dagg = 0.5690 #Total wetted area of 2 daggerboards from baseline to tip of the foil
A_rudder = 0.113 #Total immersed area of strut (1rudder)
A_elevator = 0.05 #total immersed area of elevator 1 rudder (T-part)
t_c_dagg = t_dagg/c_dagg
t_c_rudder = t_rudder/c_rudd
y0 = 1.130 #distance between baseline and tip of daggerboards
A_proj_tot = 0.26035 #projected area of 2 fwd foils fully immerged
AR_elevator = 5
form_factor_dagg = 1 + 2*(t_c_dagg) + 60*(t_c_dagg)**4
form_factor_rudder = 1 + 2*(t_c_rudder) + 60*(t_c_rudder)**4

def drag_app(Vs):
    "calculates total skin friction of immersed appendages "
    "assumes fully displacement mode"
    "assumes constant AoA for induced drag of 4°, Cl = 0.4"
    "only assumes the daggerboard induced drag for now (does not consider induced drag from rudder)"
    Re_dagg = Vs*c_dagg/nu_water
    Re_rudder = Vs*c_rudd/nu_water
    #Re_elevator = Vs*c_rudd/nu_water
    Cf_dagg = 0.075/((math.log10(Re_dagg)-2)**2)
    Cf_rudder = 0.075/((math.log10(Re_rudder)-2)**2)
    #Cf_elevator = 0.075/((math.log10(Re_elevator)-2)**2)
#    D_dagg = 2*0.5*rho_water*Vs**2*A_dagg*Cf_dagg
#    D_rudder = 2*0.5*rho_water*Vs**2*A_rudder*Cf_rudder
#    D_elevator = 2*0.5*rho_water*Vs**2*A_elevator*Cf_elevator
    Cdp_dagg = 2*(form_factor_dagg)*Cf_dagg
    Cdp_rudder = 2*(form_factor_rudder)*Cf_rudder
    Cl = 0.4
    k = 1/(2*math.pi*7.096)
    Cdi = k*Cl**2
    Cd_dagg = Cdi + Cdp_dagg
    Dp_dagg = 2*0.5*rho_water*(Vs)**2*A_dagg_up*Cd_dagg
    Dp_rudder = 2*0.5*rho_water*(Vs)**2*A_rudder*Cdp_rudder
    return Dp_dagg + Dp_rudder

def solve_lift(Fx, Zce, D_app, windage, x_crew):
    """Returns the values of lift needed for both foils to support the boat's and crews weights"""
#    x_crew = 1.1
    x_cog = 2.2
    A = np.array([[-0.08, 1.935], [1,1]])
    B = np.array([[x_crew*crew_weight + x_cog*boat_weight + Zce*Fx + 0.5*D_app - 0.5*windage] , [2707.35]])
    C = np.linalg.solve(A,B)
    L_aft = float(C[0])
    L_fwd = float(C[1])
    return (L_aft, L_fwd)

def get_lift_fwdfoil(yrh, Vs, Cl):
    """returns the lift produced by the fwd foils"""
#    A_projected = A_dagg_up*(1-yrh/y0) #distance between waterline and tip of daggerboards
    #Cl = 0.4
    L = rho_water*A_proj_tot*Cl*(Vs)**2
    return L

def get_lift_aftfoil(Vs, Cl):
    return rho_water*A_elevator*Cl*(Vs)**2


def objective(x):
    x1 = x[0] #Cl
    x2 = x[1] #alpha
    x3 = x[2] #yrh, vertical distance waterline to foil's tip
    x4 = x[3] #Boat's speed, will be fixed in the VPP loop by the bound
    x5 = x[4] #Lift required by the solve_lift function
    Re = (x4)*c_dagg/nu_water
    Cdf = 0.075/((math.log10(Re)-2)**2)
    return 0.5*rho_water*(x4)**2*(0.5576*x3-0.0315)*(2*form_factor_dagg*Cdf+ ((x1**2)/(math.pi*7.096))+(x1**2*g*c_dagg*math.exp(-g*x3/(x4)**2))/(2*(x4)**2))

def constraint1(x):
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    return 0.074*x2 + 0.116 - x1

def constraint2(x):
    x1 = x[0] #Cl
    x2 = x[1] #alpha
    x3 = x[2] #yrh, vertical distance waterline to foil's tip
    x4 = x[3] #Boat's speed Vs, will be fixed in the VPP loop by the bound
    x5 = x[4] #Lift required by the solve_lift function
    return (x5/((1/2)*rho_water*(x4)**2*(-0.4408*x3**2 + 0.8386*x3-0.1332))) -x1
    
#b1 = (-0.2, 1.2)
#b2 = (-4, 20)
#b3 = (0.3, 1.13)
#b4 = (7, 7)
#b5 = (3200, 3200)
#bnds = (b1, b2, b3, b4, b5)
con1 = {'type':'eq', 'fun':constraint1}
con2 = {'type':'eq', 'fun':constraint2}
cons = (con1, con2)

#def minimize_drag(x0, Lift_req):
#    return minimize(objective, x0, method='SLSQP', bounds=bnds, constraints=cons)

def spray_drag_dagg(Vs):
    Cdspray = 0.009 + 0.013*(t_c_dagg)
    return t_dagg*c_dagg*(Vs)**2*rho_water*Cdspray #for 2 daggerboards

def spray_drag_rudder(Vs):
    Cdspray = 0.009 + 0.013*(t_c_rudder)
    return t_rudder*c_rudd*(Vs)**2*rho_water*Cdspray #for 2 rudders

def Cdp_dagg(alpha):
#    Re_dagg = (Vs)*c_dagg/nu_water
#    Cf = 0.075/((math.log10(Re_dagg)-2)**2)
    return 0.0004*alpha**2 + 0.0008*alpha + 0.0152
    #return 2*form_factor_dagg*Cf

def Cdp_rudder(Vs):
    Re_rudder = (Vs)*c_rudd/nu_water
    Cf = 0.075/((math.log10(Re_rudder)-2)**2)
    return 2*form_factor_rudder*Cf

def Cdi_dagg(Cl):
    return Cl**2/(math.pi*7.096)

def Cdw_dagg(Vs, Cl, yrh):
    return (g*c_dagg*math.exp(-g*yrh/(2*(Vs)**2)))/(2*(Vs)**2)*Cl**2

def Cd_dagg(Cdp, Cdi, Cdw):
    return Cdp+Cdi+Cdw
   
def drag_dagg(Vs, Cl, yrh, Dspray):
    Re = (Vs)*c_dagg/nu_water
    Cdf = 0.075/((math.log10(Re)-2)**2)
    return 0.5*rho_water*(Vs)**2*(0.5576*yrh-0.0315)*(2*form_factor_dagg*Cdf+ ((Cl**2)/(math.pi*7.096))+(Cl**2*g*c_dagg*math.exp(-g*yrh/(Vs)**2))/(2*(Vs)**2))+ Dspray

def drag_rudder(Vs, Cdp, yrh, Dspray):
    A_wetted = 2*(c_rudd*(yrh-0.4)) #Wetted area for 2 rudders struts 
    return 0.5*rho_water*A_wetted*Cdp*(Vs)**2 + Dspray

def R_t(Vs, Cl_fwd, heel_angle):
    D_dagg = drag_dagg(Vs, Cd_dagg(Cdp_dagg(4),Cdi_dagg(Cl_fwd), Cdw_dagg(Vs,Cl_fwd, 1.13)),1.13, spray_drag_dagg(Vs))
    D_rudder = drag_rudder(Vs, Cdp_rudder(Vs), 1.11, spray_drag_rudder(Vs))
#    windage = get_total_windage()
    D_hull = getRh(getRf(Vs, heel_angle), getRr(Vs, heel_angle))
    return D_dagg + D_rudder + D_hull

### ---- Elevator lift & drag calculations ------ ###
    
def lift_curve_slope_rudder(Vs, yrh):
    """returns the 2d Lift curve slope Cl/Cl(inf)"""
    h = yrh/2
    k_0 = g/((Vs)**2)
    return 1/(1+(1/32)*(c_rudd/h)+math.pi*c_rudd*k_0*math.exp(-2*k_0*h))

def downwash_factor(yrh):
    h = yrh/2
    return 2.0 - 1/(math.sqrt(1+(s_rudder/(4*h))**2))

def delta_alpha_inf(Cl):
    """returns the change in AoA due to downwash"""
    return Cl/(math.pi*(s_rudder/c_rudd))

def delta_alpha(d_a_inf, downwash_factor):
    """returns the resulting change in AoA due to downwash"""
    return d_a_inf*downwash_factor

def get_AoA(Cl, lift_curve_slope, delta_alpha):
    """returns the angle of attack at a given lift coeff after downwash effects"""
    return Cl/(lift_curve_slope) + alpha_0_rudder + delta_alpha

def Cd0_elevator(alpha):
    return 0.0004*alpha**2 + 0.0008*alpha + 0.0152

def get_Cdi_elevator(Cl):
    return Cl**2/(math.pi*AR_elevator)

def get_Cdw_elevator(Cl, yrh, Vs):
    h = yrh/2
    k_0 = g/((Vs)**2)
    return 0.5*k_0*c_rudd*Cl**2*math.exp(-2*k_0*h)

def get_Cd_elevator(Cd0, Cdi, Cdw):
    return Cd0+Cdi+Cdw

def D_elevator(Cd, Vs):
    """Drag of the 2 elevators in the water"""
    return 2*0.5*rho_water*Cd*A_elevator*(Vs)**2


# In[65]:
#TWA = np.arange(25, 180, 10)
#TWS = np.arange(7, 8, 1)

def plotVPP(V_i, TWA, TWS, x_crew):
    start_time = time.clock()
    j=0
    for x in TWS:
        colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'blue', 'olive', 'cyan']
#        color_foil = ['--b', '--o', '--g', '--r', '--purple', '--brown', '--pink', '--olive', '--cyan']
        Vs1, Vb_foil, Vb_displ, theta_list, theta_foil, theta_displ, VMG = [], [], [], [], [], [], []
        for i in TWA:
            Vs = V_i
            theta_list.append(rad(i))
            (heel_angle, delta_F, Fx, Rt, c, windage, D_hull, D_dagg, D_rudder, D_elev), yrh = np.zeros(10), 1.13
            (F, R, counter) = ([],[], [])
            (Zce, flat) = (3.8, 1.0)
            loop = True
            data = {'Zce':[], 'Vt':[], 'Va': [], 'beta': [], 'Fy':[], 'HM':[], 'PM':[],
                    'Vs':[], 'deltaf':[], 'flat fact':[], 'Fx':[], 'Rt':[],
                    'Windage':[], 'D_hull':[], 'D_dagg':[], 'D_rudder':[], 'D_elevator':[], 'Immersion':[]}
            df_1 = pd.DataFrame(data)
            while loop:
                c+=1
                R.append(Rt)
                counter.append(c)
                beta = getAWA(x, i, Vs)
                Va = getAWS(x, i, Vs)
                Fx = getFx(Va, beta, i, flat)
                F.append(Fx)
                Fy = getFy(Va, beta, i, flat)
                HM = getHM1(Fy, Zce)
                PM = getPM(Fx, Zce)
                windage = get_total_windage(Va, beta)
                if HM > RM_max or PM > PM_max:
                    flat = flat - 0.01
                    Zce = getZce(Zce, flat)
                    new_row = {'Zce':Zce, 'Vt': x, 'Va': Va, 'beta': beta, 'Fy':Fy,
                               'HM':HM, 'PM':PM, 'Vs':Vs, 'deltaf': delta_F,
                               'flat fact': flat, 'Fx':Fx, 'Rt':Rt, 'Windage':windage,
                               'D_hull':D_hull, 'D_dagg': D_dagg, 'D_rudder':D_rudder, 'D_elevator':D_elev, 'Immersion':yrh}
                    df_1 = df_1.append(new_row, ignore_index = True)
                    if flat == 0.62: 
                        print("Wind too strong, sail flattened to its max")
                        break
                    else:
                        continue
                else:          
                    heel_angle = 6
                    L_req = solve_lift(Fx, Zce, drag_app(Vs), windage, x_crew)
                    Cl_fwd = L_req[1]/(0.5*A_proj_tot*rho_water*(Vs)**2)
                    #Cl_aft = L_req[0]/(0.5*2*A_elevator*rho_water*(Vs)**2)
                    if -0.2 <= Cl_fwd <= 0.8: #meaning boat can fly
                        b1 = (-0.18, 0.8) #bounds Cl_fwd
                        b2 = (-4, 9.2435) #bounds AoA
                        b3 = (0.7, 1.11) #bounds immersion, yrh
                        b4 = (Vs-delta_F, Vs+delta_F) #use boat's speed from current loop
                        b5 = (L_req[1]-0.1*L_req[1], L_req[1]+0.1*L_req[1]) #bounds Lift required for fwd foils
                        bnds = (b1, b2, b3, b4, b5) #constraints cons established beforehand
                        x0 = (Cl_fwd, 0.074*Cl_fwd+0.116, 1.1, Vs, L_req[1])
                        sol = minimize(objective, x0, method='SLSQP', bounds=bnds, constraints=cons)
                        yrh = float(sol.x[2])
                        D_dagg = sol.fun + spray_drag_dagg(Vs)
                        D_rudder = drag_rudder(Vs, Cdp_rudder(Vs), float(sol.x[2]),spray_drag_rudder(Vs))
                        D_hull = 0
                        ### -- Elevator drag calculations -- ###
                        lift_curve = lift_curve_slope_rudder(Vs, yrh)
                        d_factor = downwash_factor(yrh)
                        d_a_inf = delta_alpha_inf(Cl_fwd)
                        d_alpha = delta_alpha(d_a_inf, d_factor)
                        AoA = get_AoA(Cl_fwd, lift_curve, d_alpha)
                        Cd0 = Cd0_elevator(AoA)
                        Cdi = get_Cdi_elevator(Cl_fwd)
                        Cdw = get_Cdw_elevator(Cl_fwd, yrh, Vs)
                        Cd = get_Cd_elevator(Cd0, Cdi, Cdw)
                        Drag_elevator = D_elevator(Cd, Vs)
                        D_rudder = D_rudder + Drag_elevator
                        Rt = D_dagg + D_rudder + windage 
                    else:
                        Cl_fwd = 0.4 #Cl with best L/D ratio
                        D_dagg = drag_dagg(Vs, Cl_fwd, 1.1, spray_drag_dagg(Vs))
                        D_rudder = drag_rudder(Vs, Cdp_rudder(Vs), 1.11, spray_drag_rudder(Vs)) + 0.75*drag_rudder(Vs, Cdp_rudder(Vs), 1.13, spray_drag_rudder(Vs))
                        if abs((get_lift_fwdfoil(1.1, Vs, Cl_fwd)-L_req[1])/L_req[1]) <0.2:                
                            D_hull = get_hull_R(Vs, crew_weight, heel_angle)
        #                    D_hull = 0.7*(getRf(Vs, 5)+ getRr(Vs, 5))
                        else:
                            D_hull = get_hull_R(Vs, crew_weight, heel_angle)
        #                    D_hull = getRf(Vs, 5) + getRr(Vs, 5)
                        Rt = D_dagg + D_rudder + D_hull + windage 
        #                print('boat is not flying')
                    delta_F = (Fx - Rt)/Fx
                    if delta_F >= 0.001:
                        Vs = Vs+delta_F
                        new_row = {'Zce':Zce, 'Vt': x, 'Va': Va, 'beta': beta,
                                   'Fy':Fy, 'HM':HM, 'PM':PM, 'Vs':Vs, 'deltaf': delta_F,
                                   'flat fact': flat, 'Fx':Fx, 'Rt':Rt,
                                   'Windage':windage, 'D_hull':D_hull, 'D_dagg':D_dagg,
                                   'D_rudder':D_rudder, 'D_elevator':D_elev, 'Immersion':yrh}
                        df_1 = df_1.append(new_row, ignore_index = True)
                        continue
                    else:
                        new_row = {'Zce':Zce, 'Vt': x, 'Va': Va, 'beta': beta,
                                   'Fy':Fy, 'HM':HM, 'PM':PM, 'Vs':Vs, 'deltaf': delta_F,
                                   'flat fact': flat, 'Fx':Fx, 'Rt':Rt,
                                   'Windage':windage, 'D_hull':D_hull, 'D_dagg':D_dagg,
                                   'D_rudder':D_rudder, 'D_elevator':D_elev, 'Immersion':yrh}
                        Vs1.append(Vs/0.5144)
                        break
            cols = ['Zce', 'Vt', 'Va','beta', 'Fy', 'HM', 'PM','flat fact', 'Vs', 'Fx',
                        'Rt','Windage','D_dagg', 'D_rudder','D_hull','D_elevator', 'deltaf','Immersion']
            df_1 = df_1[cols]
            display(df_1)
#        for t in range(len(TWA)):
#            VMG.append(Vs1[t]*abs(cos(rad(TWA[t]))))
#        best_VMG_upwind = max(VMG[0:(int(len(TWA)/2))])
#        best_VMG_downwind = max(VMG[(int(len(TWA)/2)):-1])
#        print('best VMG upwind in {} knots of wind is {}knots at {}°'.format(x/0.5144, best_VMG_upwind, list(TWA)[VMG.index(best_VMG_upwind)]))
#        print('best VMG downwind in {} knots of wind is {}knots at {}°'.format(x/0.5144, best_VMG_downwind, list(TWA)[(VMG.index(best_VMG_downwind))]))
#        VMG_up, VMG_down, VMG_angle_up, VMG_angle_down=[], [], [], []
#        VMG_up.append(best_VMG_upwind)
#        VMG_down.append(best_VMG_downwind)
#        VMG_angle_up.append(list(TWA)[VMG.index(best_VMG_upwind)])
#        VMG_angle_down.append(list(TWA)[VMG.index(best_VMG_downwind)])
#        plt.figure(1)
#        plt.plot(VMG_angle_up, VMG_up, linestyle='none', marker='x', label='TWS=%s knots' %(int(x/0.5144)))
#        plt.xlabel('True wind angle (°)')
#        plt.ylabel('VMG (knots)')
#        plt.title('Upwind VMG')
#        plt.legend(loc='best')
#        plt.show()
#        plt.figure(2)
#        plt.plot(VMG_angle_down, VMG_down, linestyle='none', marker='x', label='TWS=%s knots' %(int(x/0.5144)))
#        plt.xlabel('True wind angle (°)')
#        plt.ylabel('VMG (knots)')
#        plt.title('Downwind VMG')
#        plt.legend(loc='best')
#        plt.show()
        ax = plt.subplot(projection='polar')
        ax.set_theta_offset(np.pi/2)
        ax.set_thetamin(0)
        ax.set_thetamax(180)
#        ax.set_title("Polar diagram for a 120kg crew for a range of wind speeds", fontsize = 18)
        Vs1 = interp1d(theta_list, Vs1, kind='cubic')
        plt.polar(np.linspace(rad(TWA[0]), rad(TWA[-1]), 50), Vs1(np.linspace(rad(TWA[0]), rad(TWA[-1]), 50)), ':', color = colors[j], label='TWS={}m/s - {}knots'.format(int(x), round(x/0.5144, 1)))
        theta_new1 = []
        theta_new2 = []
        Vb_displ_new1 = []
        Vb_displ_new2 = []
        if len(Vb_foil) >= 2:
            for u in range(len(theta_displ)):
                if theta_displ[u] < theta_foil[0]:
                    theta_new1.append(theta_displ[u])
                    Vb_displ_new1.append(Vb_displ[u])
                if theta_displ[u] > theta_foil[-1]:
                    theta_new2.append(theta_displ[u])
                    Vb_displ_new2.append(Vb_displ[u])
            plt.polar(theta_foil, Vb_foil, color=colors[j], linestyle='dashed')
            plt.polar(theta_new1, Vb_displ_new1, color = colors[j], label = 'TWS={}m/s - {} knots'.format(int(x), round(x/0.5144, 1)))
            plt.polar(theta_new2, Vb_displ_new2, color = colors[j])
            try:
                plt.polar([theta_new1[-1], theta_foil[0]], [Vb_displ_new1[-1],Vb_foil[0]], color=colors[j], linestyle='dashed')
                plt.polar([theta_new2[0], theta_foil[-1]], [Vb_displ_new2[0],Vb_foil[-1]], color=colors[j], linestyle='dashed') 
            except IndexError:
                continue
        else:
            print(theta_displ)
            plt.polar(theta_list, Vs1, color = colors[j], label = 'TWS={}m/s - {}knots'.format(int(x), round(x/0.5144, 1)))
        j+=1
    plt.polar(0.5, 5, label = 'foiling state', color = 'black', linestyle = 'dashed')    
    plt.legend(fontsize = 12, loc='upper left')
    print(time.clock() - start_time, "seconds")
#    plt.savefig("Polar plot 0-90° (1).png")
    return plt.show()


# In[57]:

def VPP(Vt, gamma, Vs=2, x_crew=1):
    start_time = time.clock()
    (heel_angle, delta_F, Fx, Rt, c, windage, D_hull, D_dagg, D_rudder, D_elev, yrh) = np.zeros(11)
    (Zce, flat) = (3.8, 1.0)
    loop = True
    (counter, Vs1, R, F, beta_list, Lift_list, Drag) = ([],[],[],[], [], [],[])
    data = {'Va': [], 'beta': [], 'Fy':[], 'HM':[], 'PM':[],
            'Vs':[], 'deltaf':[], 'Fx':[], 'Rt':[],
            'Windage':[], 'D_hull':[], 'D_dagg':[], 'D_rudder':[], 'D_elevator':[]}
    df_1 = pd.DataFrame(data)
    while loop:
        c+=1
        Vs1.append(Vs)
        counter.append(c)
        beta = getAWA(Vt, gamma, Vs)
        beta_list.append(beta)
        Va = getAWS(Vt, gamma, Vs)
        Fx = getFx(Va, beta, gamma, flat)
        Fy = getFy(Va, beta, gamma, flat)
        HM = getHM1(Fy, Zce)
        PM = getPM(Fx, Zce)
        R.append(Rt)
        F.append(Fx)
        windage = get_total_windage(Va, beta)
        if HM > RM_max or PM > PM_max:
            flat = flat - 0.01
            Zce = getZce(Zce, flat)
            new_row = {'Va': Va, 'beta': beta, 'Fy':Fy,
                       'HM':HM, 'PM':PM, 'Vs':Vs, 'deltaf': delta_F,
                       'Fx':Fx, 'Rt':Rt, 'Windage':windage,
                       'D_hull':D_hull, 'D_dagg': D_dagg, 'D_rudder':D_rudder}
            df_1 = df_1.append(new_row, ignore_index = True)
            if flat == 0.62: 
                print("Wind too strong, sail flattened to its max")
                break
            else:
                continue
        else:          
            heel_angle = 6
            L_req = solve_lift(Fx, Zce, drag_app(Vs), windage, x_crew)
            Cl_fwd = L_req[1]/(0.5*A_proj_tot*rho_water*(Vs)**2)
            if -0.2 <= Cl_fwd <= 0.8: #meaning boat can fly
                b1 = (-0.18, 0.8) #bounds Cl_fwd
                b2 = (-4, 9.2435) #bounds AoA
                b3 = (0.7, 1.11) #bounds immersion, yrh
                b4 = (Vs-delta_F, Vs+delta_F) #use boat's speed from current loop
                b5 = (L_req[1]-0.1*L_req[1], L_req[1]+0.1*L_req[1]) #bounds Lift required for fwd foils
                bnds = (b1, b2, b3, b4, b5) #constraints cons established beforehand
                x0 = (Cl_fwd, 0.074*Cl_fwd+0.116, 1.1, Vs, L_req[1])
                sol = minimize(objective, x0, method='SLSQP', bounds=bnds, constraints=cons)
                yrh = float(sol.x[2])
                D_dagg = sol.fun + spray_drag_dagg(Vs)
                D_rudder = drag_rudder(Vs, Cdp_rudder(Vs), float(sol.x[2]),spray_drag_rudder(Vs))
                D_hull = 0
                ### -- Elevator drag calculations -- ###
                lift_curve = lift_curve_slope_rudder(Vs, yrh)
                d_factor = downwash_factor(yrh)
                d_a_inf = delta_alpha_inf(Cl_fwd)
                d_alpha = delta_alpha(d_a_inf, d_factor)
                AoA = get_AoA(Cl_fwd, lift_curve, d_alpha)
                Cd0 = Cd0_elevator(AoA)
                Cdi = get_Cdi_elevator(Cl_fwd)
                Cdw = get_Cdw_elevator(Cl_fwd, yrh, Vs)
                Cd = get_Cd_elevator(Cd0, Cdi, Cdw)
                Drag_elevator = D_elevator(Cd, Vs)
                D_rudder = D_rudder + Drag_elevator
                Rt = D_dagg + D_rudder + windage 
            else:
                Cl_fwd = 0.4 #Cl with best L/D ratio
                D_dagg = drag_dagg(Vs, Cl_fwd, 1.1, spray_drag_dagg(Vs))
                D_rudder = drag_rudder(Vs, Cdp_rudder(Vs), 1.11, spray_drag_rudder(Vs)) + 0.75*drag_rudder(Vs, Cdp_rudder(Vs), 1.13, spray_drag_rudder(Vs))
                if abs((get_lift_fwdfoil(1.1, Vs, Cl_fwd)-L_req[1])/L_req[1]) <0.2:                
                    D_hull = get_hull_R(Vs, crew_weight, heel_angle)
#                    D_hull = 0.7*(getRf(Vs, 5)+ getRr(Vs, 5))
                else:
                    D_hull = get_hull_R(Vs, crew_weight, heel_angle)
#                    D_hull = getRf(Vs, 5) + getRr(Vs, 5)
                Rt = D_dagg + D_rudder + D_hull + windage 
#                print('boat is not flying')
            delta_F = (Fx - Rt)/Fx
#            R.append(Rt)
#            F.append(Fx)
            if delta_F >= 0.001:
                Vs = Vs+delta_F
                new_row = {'Va': Va, 'beta': beta,
                           'Fy':Fy, 'HM':HM, 'PM':PM, 'Vs':Vs, 'deltaf': delta_F,
                           'Fx':Fx, 'Rt':Rt,
                           'Windage':windage, 'D_hull':D_hull, 'D_dagg':D_dagg,
                           'D_rudder':D_rudder, 'D_elevator':D_elev}
                df_1 = df_1.append(new_row, ignore_index = True)
                continue
            else:
                new_row = {'Va': Va, 'beta': beta,
                           'Fy':Fy, 'HM':HM, 'PM':PM, 'Vs':Vs, 'deltaf': delta_F,
                           'Fx':Fx, 'Rt':Rt,'Windage':windage, 'D_hull':D_hull, 'D_dagg':D_dagg,
                           'D_rudder':D_rudder, 'D_elevator':D_elev}
                break
    cols = ['Va','beta', 'Fy', 'HM', 'PM','Vs', 'Fx',
            'Rt','Windage','D_dagg', 'D_rudder','D_hull','D_elevator', 'deltaf']
    df_1 = df_1[cols]
#    display(df_1)
#    print([D_hull, D_rudder, D_dagg, windage, Drag_elevator])
#    try:
#        print(sol)
#    except UnboundLocalError:
#        print('')
#    plt.figure(1)
#    plt.plot(beta_list, Lift_list,'-r', marker='o', markevery=10, label = 'Evolution of Lift TWA= %s°'%gamma)
#    plt.plot(beta_list, Drag, 'b', label='Evolution of Drag TWA=%s'%gamma)
#    plt.text(beta_list[0], Lift_list[0]+10, 'start')
#    plt.text(beta_list[-1], Lift_list[-1]+10, 'finish')
#    plt.plot(beta_list, F, '-y', marker='x', label='Evolution of Fx TWA=%s'%gamma)    
#    plt.plot(beta_list, [x*sin(rad(y)) - z*cos(rad(y)) for x,y,z in zip(Lift_list, beta_list, Drag)], '-y', label='Fx formula output')
#    plt.plot(counter, R, '-r', label = 'Overall Resistance')
#    plt.legend(loc='best', fontsize = 17)
#    plt.xlabel('AWA', fontsize = 17)
#    plt.ylabel('Forces, N', fontsize = 17)
#    plt.title('Lift, Drag and driving Forces vs. AWA', fontsize = 21)
#    plt.show()
#    print(R)
#    print(F)
#    print(D_hull)
#    print('xcrew = {}:final boat speed for TWA = {}° and TWS = {} knots is : {} m/s ({} knots)'.format(x_crew, gamma, round(Vt/0.5144, 1), round(Vs, 3), round(Vs/0.5144, 3)))
#    print(Vs)
    print(time.clock() - start_time, "seconds")
    print('Final boat speed is:{}m/s ({} knots)'.format(round(Vs, 3), round(Vs/0.5144,3)))
    return (Vs, df_1)



# In[58]:
        
def f(x_crew, TWA_range, TWS_range):
    start_time = time.clock()
    """Polar diagram with optimised crew position and information on foiling state"""
    colors = ['w', 'b', 'g', 'r', 'cyan', 'm', 'y', 'k']
    j=0
    size_df=0
    for i in TWS_range:
        VMG, theta_list, V, Vb_foil, Vb_displ, theta_foil, theta_displ, theta_new1, theta_new2,Vb_displ_new1, Vb_displ_new2 =[],[],[],[],[],[],[],[],[],[],[]
        j+=1
        for u in TWA_range:
            theta_list.append(rad(u))
            V1=[]
            for x in x_crew:
                a = VPP(i, u, 2, x)
                size_df+=a[1].size
                V1.append(round(a[0], 3))
            if int(VPP(i, u, 2, x_crew[V1.index(max(V1))])[1]['D_hull'][-1:]) == 0:
                Vb_foil.append(max(V1))
                theta_foil.append(rad(u))
            else:
                Vb_displ.append(max(V1))
                theta_displ.append(rad(u))
            V.append(max(V1))
        if len(Vb_foil) >= 2:
            for u in range(len(theta_displ)):
                if theta_displ[u] < theta_foil[0]:
                    theta_new1.append(theta_displ[u])
                    Vb_displ_new1.append(Vb_displ[u])
                if theta_displ[u] > theta_foil[-1]:
                    theta_new2.append(theta_displ[u])
                    Vb_displ_new2.append(Vb_displ[u])
            ax = plt.subplot(projection='polar')
            ax.set_theta_offset(np.pi/2)
            ax.set_thetamin(TWA_range[0]-5)
            ax.set_thetamax(TWA_range[-1]+5)
            ax.text(rad(20), 4, 'Boat speed (m/s)', fontsize=12, rotation=-60)           
            x1new=np.linspace(theta_new1[0], theta_new1[-1], 50)
            x2new=np.linspace(theta_foil[0], theta_foil[-1], 50)
            x3new=np.linspace(theta_new2[0], theta_new2[-1], 50)
            try:
                f1 = interp1d(theta_new1, Vb_displ_new1, kind='quadratic', fill_value='interpolate')
            except ValueError:
                theta_new1.append(theta_new1[-1]+0.001)
                Vb_displ_new1.append(Vb_displ_new1[-1]+0.001)
                f1 = interp1d(theta_new1, Vb_displ_new1, kind='quadratic', fill_value='interpolate')
            try:
                f2 = interp1d(theta_foil, Vb_foil, kind='quadratic', fill_value='interpolate')
            except ValueError:
                theta_foil.append(theta_foil[-1]+0.001)
                Vb_foil.append(Vb_foil[-1]+0.001)
                f2 = interp1d(theta_foil, Vb_foil, kind='quadratic', fill_value='interpolate')
            try:                
                f3 = interp1d(theta_new2, Vb_displ_new2, kind='quadratic', fill_value='interpolate')
            except ValueError:
                theta_new2.append(theta_new2[-1]+0.001)
                Vb_displ_new2.append(Vb_displ_new1[-1]+0.001)   
                f3 = interp1d(theta_new2, Vb_displ_new2, kind='quadratic', fill_value='interpolate')
            plt.plot(x1new, f1(x1new),'-', color=colors[j], label='TWS={}m/s'.format(i))
            plt.plot(x2new, f2(x2new), linestyle='dashed', color=colors[j])
            plt.plot(x3new, f3(x3new), '-', color=colors[j])
            try:
                plt.polar([theta_new1[-1], theta_foil[0]], [Vb_displ_new1[-1],Vb_foil[0]], color=colors[j], linestyle='dashed')
                plt.polar([theta_new2[0], theta_foil[-1]], [Vb_displ_new2[0],Vb_foil[-1]], color=colors[j], linestyle='dashed') 
            except IndexError:
                continue          
        else:
            f0 = interp1d(theta_list, V, kind='quadratic')
            x0new = np.linspace(theta_list[0], theta_list[-1], 100)            
            plt.plot(x0new, f0(x0new), '-', color = colors[j], label='TWS={}m/s'.format(i))
#        for t in range(len(TWA_range)):
#            VMG.append(V[t]*abs(cos(rad(TWA_range[t]))))
#        best_VMG_upwind = max(VMG[0:int(np.where(TWA_range==85)[0])])
#        best_VMG_downwind = max(VMG[int(np.where(TWA_range==105)[0]):-1])
#        print('best VMG upwind in {} knots of wind is {}knots at {}°'.format(i/0.5144, best_VMG_upwind, list(TWA_range)[VMG.index(best_VMG_upwind)]))
#        print('best VMG downwind in {} knots of wind is {}knots at {}°'.format(i/0.5144, best_VMG_downwind, list(TWA_range)[(VMG.index(best_VMG_downwind))]))
    print(time.clock() - start_time, "seconds")
    plt.polar(0.5, 5, label = 'foiling state', color = 'black', linestyle = 'dashed')    
    plt.legend(loc='upper left', fontsize = 12)
    return plt.show()



def VPP2(Vt, TWA_range, Vs, xcrew):
    colors = ['w', 'b', 'g', 'r', 'cyan', 'm', 'y', 'k']
    j=0    
#    V_max=[]
    for x in xcrew:
        V=[]
        j+=1
        for u in TWA_range:
            V.append(VPP(Vt, u, Vs, x)[0])
#        V_max.append(max(V))    
        plt.plot(TWA_range, V, marker='x', color=colors[j], linestyle='-', label='x_crew=%s m'%x)
#    plt.plot(TWA_range, V_max, '-', label='Final output')
    plt.legend(loc='best')
    plt.xlabel('TWA (°)')
    plt.ylabel('Boat speed')
    return plt.show()

def f2(Vt, TWA_range, Vs, xcrew):
    start_time = time.clock()
    colors = ['w', 'b', 'g', 'r', 'cyan', 'm', 'y', 'k']
    j=0
    for x in xcrew:
        VMG, theta_list, V, Vb_foil, Vb_displ, theta_foil, theta_displ, theta_new1, theta_new2,Vb_displ_new1, Vb_displ_new2 =[],[],[],[],[],[],[],[],[],[],[]
        j+=1
        V1=[]
        for u in TWA_range:
            theta_list.append(rad(u))
            V1.append(round(VPP(Vt, u, Vs, x)[0], 3))
#            if int(VPP(Vt, u, 2, xcrew[V1.index(max(V1))])[1]['D_hull'][-1:]) == 0:
#                Vb_foil.append(max(V1))
#                theta_foil.append(rad(u))
#            else:
#                Vb_displ.append(max(V1))
#                theta_displ.append(rad(u))
#        if len(Vb_foil) >= 2:
#            for u in range(len(theta_displ)):
#                if theta_displ[u] < theta_foil[0]:
#                    theta_new1.append(theta_displ[u])
#                    Vb_displ_new1.append(Vb_displ[u])
#                if theta_displ[u] > theta_foil[-1]:
#                    theta_new2.append(theta_displ[u])
#                    Vb_displ_new2.append(Vb_displ[u])
        ax = plt.subplot(projection='polar')
        ax.set_theta_offset(np.pi/2)
        ax.set_thetamin(TWA_range[0]-5)
        ax.set_thetamax(TWA_range[-1]+5)
        ax.text(rad(20), 3, 'Boat speed (m/s)', fontsize=12, rotation=-60)
#            x1new=np.linspace(theta_new1[0], theta_new1[-1], 50)
#            x2new=np.linspace(theta_foil[0], theta_foil[-1], 50)
#            x3new=np.linspace(theta_new2[0], theta_new2[-1], 50)
#            try:
#                f1 = interp1d(theta_new1, Vb_displ_new1, kind='quadratic', fill_value='interpolate')
#            except ValueError:
#                theta_new1.append(theta_new1[-1]+0.001)
#                Vb_displ_new1.append(Vb_displ_new1[-1]+0.001)
#                f1 = interp1d(theta_new1, Vb_displ_new1, kind='quadratic', fill_value='interpolate')
#            try:
#                f2 = interp1d(theta_foil, Vb_foil, kind='quadratic', fill_value='interpolate')
#            except ValueError:
#                theta_foil.append(theta_foil[-1]+0.001)
#                Vb_foil.append(Vb_foil[-1]+0.001)
#                f2 = interp1d(theta_foil, Vb_foil, kind='quadratic', fill_value='interpolate')
#            try:                
#                f3 = interp1d(theta_new2, Vb_displ_new2, kind='quadratic', fill_value='interpolate')
#            except ValueError:
#                theta_new2.append(theta_new2[-1]+0.001)
#                Vb_displ_new2.append(Vb_displ_new1[-1]+0.001)   
#                f3 = interp1d(theta_new2, Vb_displ_new2, kind='quadratic', fill_value='interpolate')
#            plt.plot(x1new, f1(x1new),'-', color=colors[j], label='x_crew=%s m'%x)
#            plt.plot(x2new, f2(x2new), linestyle='dashed', color=colors[j])
#            plt.plot(x3new, f3(x3new), '-', color=colors[j])
#            try:
#                plt.polar([theta_new1[-1], theta_foil[0]], [Vb_displ_new1[-1],Vb_foil[0]], color=colors[j], linestyle='dashed')
#                plt.polar([theta_new2[0], theta_foil[-1]], [Vb_displ_new2[0],Vb_foil[-1]], color=colors[j], linestyle='dashed') 
#            except IndexError:
#                continue          
#        else:
#            f0 = interp1d(theta_list, V, kind='quadratic')
#            x0new = np.linspace(theta_list[0], theta_list[-1], 100)            
        plt.plot(theta_list, V1, '-', marker='x', color = colors[j], label='x_crew=%s m'%x)
#        for t in range(len(TWA_range)):
#            VMG.append(V[t]*abs(cos(rad(TWA_range[t]))))
#        best_VMG_upwind = max(VMG[0:int(np.where(TWA_range==85)[0])])
#        best_VMG_downwind = max(VMG[int(np.where(TWA_range==105)[0]):-1])
#        print('best VMG upwind in {} knots of wind is {}knots at {}°'.format(i/0.5144, best_VMG_upwind, list(TWA_range)[VMG.index(best_VMG_upwind)]))
#        print('best VMG downwind in {} knots of wind is {}knots at {}°'.format(i/0.5144, best_VMG_downwind, list(TWA_range)[(VMG.index(best_VMG_downwind))]))
    print(time.clock() - start_time, "seconds")
#    plt.polar(0.5, 5, label = 'foiling state', color = 'black', linestyle = 'dashed')    
    plt.legend(loc='upper left', fontsize = 12)
    plt.show()
    return print('The list of optimum speeds for a range of angles between 25° and 70° is:{}'.format(V))
    

# =============================================================================
# 
#                 theta_foil.append(theta_foil[-1]+0.0001), theta_foil.append(theta_foil[-1]+0.0001)
#                 Vb_foil.append(Vb_foil[-1]+0.001), Vb_foil.append(Vb_foil[-1]+0.001)
#                 theta_new1.append(theta_new1[-1]+0.0001)
#                 Vb_displ_new1.append(Vb_displ_new1[-1]+0.001)
#                 theta_new2.append(theta_new2[-1]+0.0001)
#                 Vb_displ_new2.append(Vb_displ_new2[-1]+0.001)            
#                 V_foil = interp1d(theta_foil, Vb_foil, kind='cubic')
#                 V_displ_1 = interp1d(theta_new1, Vb_displ_new1, kind='cubic')
#                 V_displ_2 =interp1d(theta_new2, Vb_displ_new2, kind='cubic')
#     #        plt.polar(theta_new1, Vb_displ_new1, color = colors[j], label = 'TWS={}m/s - {} knots'.format(int(i), round(i/0.5144, 1)))
#     #        plt.polar(theta_foil, Vb_foil, color = colors[j], linestyle='dashed')
#     #        plt.polar(theta_new2, Vb_displ_new2, color = colors[j])
#             plt.polar(np.linspace(rad(theta_foil[0]),rad(theta_foil[-1]),50),V_foil(np.linspace(rad(theta_foil[0]), rad(theta_foil[-1]), 50)), color=colors[j], linestyle='dashed', marker='x', markevery=10)
#             plt.polar(np.linspace(rad(theta_new1[0]),rad(theta_new1[-1]),50),V_displ_1(np.linspace(rad(theta_new1[0]), rad(theta_new1[-1]), 50)), color=colors[j], linestyle='-', markevery=10, label = 'TWS={}m/s - {} knots'.format(int(i), round(i/0.5144, 1)))
#             plt.polar(np.linspace(rad(theta_new2[0]),rad(theta_new2[-1]),50),V_displ_2(np.linspace(rad(theta_new2[0]), rad(theta_new2[-1]), 50)), color=colors[j], linestyle='-', markevery=10)
#             try:
#                 plt.polar([theta_new1[-1], theta_foil[0]], [Vb_displ_new1[-1],Vb_foil[0]], color=colors[j], linestyle='dashed')
#                 plt.polar([theta_new2[0], theta_foil[-1]], [Vb_displ_new2[0],Vb_foil[-1]], color=colors[j], linestyle='dashed') 
#             except IndexError:
#                 continue
#     #        plt.polar(np.linspace(rad(TWA_range[0]), rad(TWA_range[-1]), 50), V(np.linspace(rad(TWA_range[0]), rad(TWA_range[-1]), 50)), '-', color=colors[j])
#             plt.polar(0.5, 5, label = 'foiling state', color = 'black', linestyle = 'dashed')
#         else:
#             plt.polar(np.linspace(rad(TWA_range[0]), rad(TWA_range[-1]), 50), V(np.linspace(rad(TWA_range[0]), rad(TWA_range[-1]), 50)), '-', color=colors[j])
# 
# =============================================================================

"""
                        Parabolic Solar Dish Stirling System
                              By Pablo Jimenez Zabalaga
"""
import math
import numpy as np
import matplotlib.pyplot as plt

# Input data
# Time (h)
T = np.arange(0,24)
# Direct Normal Irradiation (W/m²)
DNI = [0, 0, 0, 0, 0, 0, 0, 100, 200, 350, 500, 650, 900, 850, 650, 450, 300, 200, 100, 0, 0, 0, 0, 0] # (W/m^2) 
# Wind speed (m/s)
V_air = [1, 1, 1, 2, 2, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] # (m/s)
# Ambient temperature (ºC)
T_amb = [15, 14, 13, 12, 13, 14, 15, 17, 19, 21, 22, 23, 25, 26, 27, 28, 27, 25, 23, 21, 19, 17, 16, 15]
# Area of the concentrator (m²)
A_con = 56.7
# Optical efficiency of concentrator (%)
eta_op = 94
# Interception factor (%)
phi_fi = 90
# Air density (kg/m³)
rho_air = 1.225
# Diameter of receiver (m)
D_rec = 0.30
# Area of receiver (m²)
A_rec = D_rec**2*(math.pi/4)
# Aperture diameter of cavity (m)
D_ap = 0.19
# Aperture area of cavity (m²)
A_ap = D_ap**2*(math.pi/4)
# Concentration factor
C = A_con/A_ap
# Dynamic viscosity of air (kg/(m*s))
mu_air = 1.802e-5
# Thermal conductivity of air (W/K*m)
lambda_air = 0.025
# Receiver absorptivity coefficient
alpha_rec = 0.95
# Stefan-Boltzmann constant (W/(m²*K⁴))
sigma = 5.67e-8
# Absorber emissivity
epsilon_ab = 0.9
# Deterioration factor
F_d = 90
# Cooling factor
R_ref = 50
# Attenuation constant
K_att = (eta_op/100)*(phi_fi/100)*(F_d/100)*(R_ref/100)
# Thickness of the ceramic wall (m)
e_cer = 0.20
# Thermal conductivity of ceramic (W/(K*m))
lambda_cer = 10
# Gravity (m/s²)
g = 9.81
# Coefficient of thermal expansion of air (1/K)
beta = 0.0034
# Depth of the cavity (m)
L_rec = 0.12
# Kinematic viscosity of air (m²/s)
upsilon_air = 1.48e-5
# Parameter of Equation
S = 1.12-0.982*(D_ap/L_rec)
# Inclination angle (º)
theta_inc = 40
# function
f_theta_inc = 0.1634 +0.7498*math.sin(math.radians(theta_inc))-0.5026*math.sin(math.radians(2*theta_inc))+0.3278*math.sin(math.radians(3*theta_inc))
# Absorption coefficient of the receiving cavity
alpha_eff = alpha_rec/(alpha_rec+(1-alpha_rec)*(A_ap/A_rec))
# Compression space temperature (°C)
T_c = 25
# Expansion space temperature (°C)
T_e = 850
# Temperature correction term
T_correct = 1-math.sqrt((T_c+273.15)/(T_e+273.15))
# Mean engine pressure (Pa)
P_mean = 15000000
# Swept volume of the engine (m³)
V_sw = 0.00016
# Engine frequency (Hz)
f = 25
# Power received by the concentrator (W)
Q̇_s = []
# Power received by the receiver (W)
Q̇_rec = []
# Reynolds number
Re = []
# Natural convection heat transfer coefficient of receiver (W/(m²*K))
h_ext = []
# Ambient temperature (K)
T_amb_K = []
# Absorber temperature (K)
T_ab_K = []
# Conduction losses (W)
Q̇_cond = []
# Grashof number
Gr = []
# Nusselt number
Nu = []
# Natural convection heat transfer coefficient (W/(m²*K))
h_nat = []
# Forced convection heat transfer coefficient (W/(m²*K))
h_forced = []
# Total convection heat transfer coefficient (W/(m²*K))
h_total = []
# Convection losses (W)
Q̇_conv = []
# Radiation losses (W)
Q̇_rad = []
# Reflection losses (W)
Q̇_ref = []
# Power intercepted by receiver (W)
Q̇_u = []
# Beale curve
Beale_curve = []
# Gross power of Stirling engine
P_gross = []


# Output data
for t in T:
    # Power received by the concentrator (W)
    Q̇_s.append(A_con*DNI[t])
    
    # Power received by the receiver (W)
    Q̇_rec.append(Q̇_s[t]*(eta_op/100)*(phi_fi/100))
    
    # Conduction losses (W)
    # Reynolds number
    Re.append((rho_air*V_air[t]*D_rec)/mu_air)
    
    # Natural convection heat transfer coefficient of receiver (W/(m²*K))
    h_ext.append((0.148*lambda_air*Re[t]**0.633)/D_rec)
    
    # Ambient temperature (K)
    T_amb_K.append(T_amb[t]+273.15)
    
    # Absorber temperature (K)
    if(((DNI[t]*C*alpha_rec)/(sigma*epsilon_ab))**(1/4)*K_att==0):
        T_ab_K.append(T_amb_K[t])
    else:
        T_ab_K.append(((DNI[t]*C*alpha_rec)/(sigma*epsilon_ab))**(1/4)*K_att)
    
    # Conduction losses (W)
    Q̇_cond.append((T_ab_K[t]-T_amb_K[t])/((e_cer/(lambda_cer*A_rec))+(1/(h_ext[t]*A_rec))))
    
    # Grashof number
    Gr.append((g*beta*(L_rec**3)*(T_ab_K[t]-T_amb_K[t]))/(upsilon_air**2))
    
    # Nusselt number
    Nu.append(0.088*(Gr[t]**(1/3))*((math.cos(math.radians(theta_inc)))**2.47)*((T_ab_K[t]/T_amb_K[t])**0.18)*((D_ap/L_rec)**S))
    
    # Natural convection heat transfer coefficient (W/(m²*K))
    h_nat.append((Nu[t]*rho_air)/L_rec)
    
    # Forced convection heat transfer coefficient (W/(m²*K))
    h_forced.append(f_theta_inc*(V_air[t]**1.401))
    
    # Total convection heat transfer coefficient (W/(m²*K))
    h_total.append(h_nat[t]+h_forced[t])
    
    # Convection losses (W)
    Q̇_conv.append(h_total[t]*A_ap*(T_ab_K[t]-T_amb_K[t]))
    
    # Radiation losses (W)
    Q̇_rad.append(epsilon_ab*sigma*A_ap*((T_ab_K[t]**4)-(T_amb_K[t]**4)))
    
    # Reflection losses (W)
    Q̇_ref.append((1-alpha_eff)*Q̇_rec[t])
    
    # Power intercepted by receiver (W)
    Q̇_u.append(Q̇_rec[t]-Q̇_ref[t]-Q̇_cond[t]-Q̇_conv[t]-Q̇_rad[t])
    
    # Beale curve
    Beale_curve.append(0.0391+0.0000221*Q̇_u[t]-3.55e-10*(Q̇_u[t]**2))
    
    # Gross power of Stirling engine
    P_gross.append(Beale_curve[t]*P_mean*V_sw*f*T_correct)
    
# Graph
P_gross_fig, ax1 = plt.subplots(linewidth=6, edgecolor="k")
ax1.set_xlabel('Time (h)',size='14' ,weight='normal') # x-axis parameter labeling
ax1.set_ylabel('P$_{gross}$ (W)',size='14',weight='normal')
ax1.plot(T,P_gross,linewidth=1,color='b',label='P$_{gross}$ (W)')
ax1.tick_params(axis='y')
ax1.grid(color='b',linestyle=':',linewidth='0.5') # grid on
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set_ylabel('DNI (W/m$^2$)',size='14',weight='normal')  # we already handled the x-label with ax1
ax2.plot(T,DNI,linewidth=1,color="g",label='DNI (W/m$^2$)')
ax2.tick_params(axis='y')
P_gross_fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)
P_gross_fig.tight_layout()  # otherwise the right y-label is slightly clipped
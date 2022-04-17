import numpy as np
import pandas as pd 
from pandas import DataFrame
from scipy.interpolate import interp1d
import scipy.integrate as integrate
import math

############################ Reactor Spectrum Data #########################################

# From P. Huber:https://arxiv.org/abs/1106.0687
#(²³⁹Pu, ²⁴¹Pu, ²³⁵U) - E>2MeV
data_Pu239 = pd.read_csv("Reactor/Pu239-anti-neutrino-flux-250keV.csv")
data_Pu241 = pd.read_csv("Reactor/Pu241-anti-neutrino-flux-250keV.csv")
data_U235 = pd.read_csv("Reactor/U235-anti-neutrino-flux-250keV.csv")

#Energy bins
E_P_Huber=data_Pu239.loc[:,'energy [MeV]'].values
#Pu239
Nu_Pu239_P_Huber=data_Pu239.loc[:,'neutrinos[MeV^-1 fission^-1] '].values
Nu_Pu239_P_Huber_error=[]
Nu_Pu239_P_Huber_error.append(np.absolute(data_Pu239.loc[:,'neg total error'].values)*Nu_Pu239_P_Huber)
Nu_Pu239_P_Huber_error.append(np.absolute(data_Pu239.loc[:,'pos total error'].values)*Nu_Pu239_P_Huber)
#Pu241
Nu_Pu241_P_Huber=data_Pu241.loc[:,'neutrinos[MeV^-1 fission^-1] '].values
Nu_Pu241_P_Huber_error=[]
Nu_Pu241_P_Huber_error.append(np.absolute(data_Pu241.loc[:,'neg total error'].values)*Nu_Pu241_P_Huber)
Nu_Pu241_P_Huber_error.append(np.absolute(data_Pu241.loc[:,'pos total error'].values)*Nu_Pu241_P_Huber)
#U235
Nu_U235_P_Huber=data_U235.loc[:,'neutrinos[MeV^-1 fission^-1] '] .values
Nu_U235_P_Huber_error=[]
Nu_U235_P_Huber_error.append(np.absolute(data_U235.loc[:,'neg total error'].values)*Nu_U235_P_Huber)
Nu_U235_P_Huber_error.append(np.absolute(data_U235.loc[:,'pos total error'].values)*Nu_U235_P_Huber)

#From Th. A.Muller, et.al: https://arxiv.org/abs/1101.2663v3
#(²³⁸U) - E>2MeV
data_U238 = pd.read_csv("Reactor/U238-anti-neutrino-flux-Muller.csv")
#Energy bins
E_Muller=data_U238.loc[:,'Kinetic_E(MeV)'].values
#U238
Nu_U238_Muller=data_U238.loc[:,'N_nu_bar_e(/fission/MeV)-450d'].values

#From P. Vogel
#(²³⁹Pu, ²⁴¹Pu, ²³⁵U, ²³⁸U) - E<2MeV
E_less_2MeV= [2,1.5,1,0.75,0.5,0.25,0.125,6.25e-2, 3.12e-2, 1.563e-2, 7.813e-3,0]
U_235_less_2MeV=[1.26,1.69,2.41,2.66,2.66,2.16,1.98,0.61,0.35,0.092,0.024,0]
Pu_239_less_2MeV=[1.08, 1.48, 2.32, 2.58, 2.63, 2.08, 1.99, 0.64, 2.13, 0.56, 0.14,0]
U_238_less_2MeV=[1.5, 1.97, 2.75, 2.96, 2.91, 2.18, 2.02, 0.65, 1.32, 0.35, 0.089,0]
Pu_241_less_2MeV=[1.32, 1.75, 2.63, 2.9, 2.82, 2.14, 1.85, 0.59, 3, 0.79, 0.2,0]

################################## Interpolation ###############################################
U_235_more_2MeV_interpol = interp1d(E_P_Huber, Nu_U235_P_Huber, kind='cubic') #interpolation
Pu_239_more_2MeV_interpol = interp1d(E_P_Huber, Nu_Pu239_P_Huber, kind='cubic') #interpolation
U_238_more_2MeV_interpol = interp1d(E_P_Huber, Nu_U238_Muller, kind='cubic') #interpolation
Pu_241_more_2MeV_interpol = interp1d(E_P_Huber, Nu_Pu241_P_Huber, kind='cubic') #interpolation

U_235_less_2MeV_interpol = interp1d(E_less_2MeV, U_235_less_2MeV, kind='cubic') #interpolation
Pu_239_less_2MeV_interpol = interp1d(E_less_2MeV, Pu_239_less_2MeV, kind='cubic') #interpolation
U_238_less_2MeV_interpol = interp1d(E_less_2MeV, U_238_less_2MeV, kind='cubic') #interpolation
Pu_241_less_2MeV_interpol = interp1d(E_less_2MeV, Pu_241_less_2MeV, kind='cubic') #interpolation


def dN_dE_Muller_Huber_Vogel(E,element):
  if element=='U235':
    if E<=2.0:
      return U_235_less_2MeV_interpol(E)
    else:
     return U_235_more_2MeV_interpol(E)

  if element=='U238':
    if E<=2.0:
      return U_238_less_2MeV_interpol(E)
    else:
     return U_238_more_2MeV_interpol(E)
     
  if element=='Pu239':
    if E<=2.0:
      return Pu_239_less_2MeV_interpol(E)
    else:
     return Pu_239_more_2MeV_interpol(E)

  if element=='Pu241':
    if E<=2.0:
      return Pu_241_less_2MeV_interpol(E)
    else:
     return Pu_241_more_2MeV_interpol(E)

dN_dE_Muller_Huber_Vogel_vec=np.vectorize(dN_dE_Muller_Huber_Vogel)


def dN_dE_Muller_Huber_Vogel_tot(E):
  U_235=dN_dE_Muller_Huber_Vogel_vec(E,'U235')
  U_238=dN_dE_Muller_Huber_Vogel_vec(E,'U238')
  Pu_239=dN_dE_Muller_Huber_Vogel_vec(E,'Pu239')
  Pu_241=dN_dE_Muller_Huber_Vogel_vec(E,'Pu241')
  dN_dE=(0.56*U_235)+(0.08*U_238)+(0.30*Pu_239)+(0.06*Pu_241)
  return dN_dE

dN_dE_Muller_Huber_Vogel_tot_vec=np.vectorize(dN_dE_Muller_Huber_Vogel_tot)



############################# Flux at the Angra Detector ################################
Norm_Muller_Huber_Vogel=integrate.quad(lambda E: dN_dE_Muller_Huber_Vogel_tot(E), 0, 8)
def Flux_Muller_Huber_Vogel(Enu): #E[MeV]
  D=30*10**2 #cm
  N_rate= 8.7*10**20# Neutrino rate of production (s⁻¹) from CONNIE https://arxiv.org/abs/1608.01565
  Flux=(1/(4*math.pi*D**2))*N_rate*(dN_dE_Muller_Huber_Vogel_tot(Enu)/Norm_Muller_Huber_Vogel[0])#cm⁻².s⁻¹.MeV⁻¹
  return Flux 
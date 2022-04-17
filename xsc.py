import numpy as np
import scipy.integrate as integrate
import math

#Constants
G_F=1.16637*(10**(-11)) #MeV-²
sin2_W=0.23149
m_nucleon= 939.56536 #Nutron mass - MeV

#Units change
MeV_to_inv_m = 5.06773*(10**12) # From MeV to m^-1
MeV_to_inv_cm=5.06773*(10**10) #MeV -> cm^-1
inv_MeV2_to_cm2= MeV_to_inv_cm**(-2) #MeV^-2 -> cm^2 
MeV_to_inv_fm = MeV_to_inv_m*(10**(-15)) # From MeV to fm^-1

#Nuclear Form Factor
def form_f2(T,Z,N): 
  A=Z+N
  M=A*m_nucleon
  q=np.sqrt(2*M*T) #MeV
  q=q*MeV_to_inv_fm #MeV -> fm^-1
  a=0.52 #fm
  s=0.9 #fm
  c=(1.23*(A**(1/3))-0.6) #fm
  r=np.sqrt(c**2+(7/3)*(np.pi**2)*(a**2)-5*(s**2))#fm
  qr=q*r
  qs=q*s
  #print(qr)
  form=3*((np.sin(qr)-qr*np.cos(qr))/(qr**3))*np.exp(-1*((qs)**2)/2)
  #print(form**2)
  return form**2

#CEvNS Differential Cross-section
def CS_coherent(T,Enu,Z,N):
  A=Z+N
  M=A*m_nucleon
  Enu_min=np.sqrt(M*T/2)
  cs=0
  if Enu>Enu_min:
    cs=(G_F**2)*(1/(4*math.pi))*((Z*(4*sin2_W-1)+N)**2)*M*(1-(M*T/(2*(Enu**2))))*form_f2(T,Z,N) # MeV-2/MeV
    cs=cs*inv_MeV2_to_cm2 #MeV^-2 /MeV -> cm^2/MeV
  return cs 
CS_coherent_vec=np.vectorize(CS_coherent)

#Differential Cross-section for each mineral
def CS_coherent_mineral(T,Enu,mineral):
    if mineral =='SiO2':
        CS_O=CS_coherent_vec(T,Enu,8,8)
        CS_Si=CS_coherent_vec(T,Enu,14,14)
        CS=CS_Si+2*CS_O
        return CS
    if mineral =='MgSO4_7_H2O':
        CS_O=CS_coherent_vec(T,Enu,8,8)
        CS_S=CS_coherent_vec(T,Enu,16,16)
        CS_Mg=CS_coherent_vec(T,Enu,12,12)
        CS_H=CS_coherent_vec(T,Enu,1,1)
        CS= CS_Mg+CS_S+4*CS_O+7*(2*CS_H+CS_O)
        return CS
    if mineral =='NaCl':
        CS_Na=CS_coherent_vec(T,Enu,11,12)
        CS_Cl=CS_coherent_vec(T,Enu,17,18)
        CS= CS_Na+CS_Cl
        return CS

#Total Cross-section integrated in the recoil energy T
def CS_coherent_tot(Enu,Z,N): 
  A=Z+N
  M=A*m_nucleon
  T_max=2*(Enu**2)/M
  cs_tot=integrate.quad(lambda T: CS_coherent_vec(T,Enu,Z,N), 0, 2*T_max)[0]
  return cs_tot  #cm^2
CS_coherent_tot_vec=np.vectorize(CS_coherent_tot)

#Total Cross-section for each mineral
def CS_coherent_mineral_tot(Enu,mineral):
    if mineral =='SiO2':
        CS_O=CS_coherent_tot_vec(Enu,8,8)
        CS_Si=CS_coherent_tot_vec(Enu,14,14)
        CS=CS_Si+2*CS_O
        return CS
    if mineral =='MgSO4_7_H2O':
        CS_O=CS_coherent_tot_vec(Enu,8,8)
        CS_S=CS_coherent_tot_vec(Enu,16,16)
        CS_Mg=CS_coherent_tot_vec(Enu,12,12)
        CS_H=CS_coherent_tot_vec(Enu,1,1)
        CS= CS_Mg+CS_S+4*CS_O+7*(2*CS_H+CS_O)
        return CS
    if mineral =='NaCl':
        CS_Na=CS_coherent_tot_vec(Enu,11,12)
        CS_Cl=CS_coherent_tot_vec(Enu,17,18)
        CS= CS_Na+CS_Cl
        return CS

# Approximation for comparison
# def CS_coherent_tot_aprox(Enu,Z,N):
#   return 4.22*(10**(-45))*(N**2)*(Enu**2) #cm^2
from xsc import*
from Reactor_Flux import*
import scipy.integrate as integrate
from WIMpy import DMUtils as DMU

#Constants
SiO2_molar_mass=60.07 #g/mol
MgSO4_7_H2O_molar_mass= 246.4746 #g/mol
NaCl_molar_mass=58.433 #g/mol


#################################### Reactor Neutrino Signal #########################################

### Per Mineral ###
def Event_rate(Enu,T,mineral):
  R_per_target=CS_coherent_mineral(T,Enu,mineral)*Flux_Muller_Huber_Vogel(Enu)  #s⁻¹MeV⁻²

  if mineral =='SiO2':
    n_targets_density= (6.02*10**23)/SiO2_molar_mass #SiO2 targets/g
  elif mineral =='MgSO4_7_H2O':
    n_targets_density= (6.02*10**23)/MgSO4_7_H2O_molar_mass #MgSO4 targets/g
  elif mineral =='NaCl':
    n_targets_density= (6.02*10**23)/NaCl_molar_mass #MgSO4 targets/g
  else:
    print("Invalid Mineral!")
    return 0
  
  R=n_targets_density*R_per_target#s⁻¹MeV⁻²g⁻¹
  return R
     

def Recoil_spectrum(T,mineral):
  Recoil=integrate.quad(lambda Enu: Event_rate(Enu,T,mineral),0,8,epsabs=1.49e-20) #s⁻¹MeV⁻¹g⁻¹
  #print(Recoil)
  return Recoil #s⁻¹MeV⁻¹g⁻¹

### Per Atom ###
def Event_rate_atom(Enu,T,atom):
  if atom =='Si':
    R_per_target=CS_coherent(T,Enu,14,14)*Flux_Muller_Huber_Vogel(Enu)  #s⁻¹MeV⁻²
    R=(6.02*10**23)*R_per_target#s⁻¹MeV⁻²mol⁻¹
    return R

  elif atom =='O':
    R_per_target=CS_coherent(T,Enu,8,8)*Flux_Muller_Huber_Vogel(Enu)  #s⁻¹MeV⁻²
    R=(6.02*10**23)*R_per_target#s⁻¹MeV⁻²mol⁻¹
    return R

  elif atom =='S':
    R_per_target=CS_coherent(T,Enu,16,16)*Flux_Muller_Huber_Vogel(Enu)  #s⁻¹MeV⁻²
    R=(6.02*10**23)*R_per_target#s⁻¹MeV⁻²mol⁻¹
    return R

  elif atom =='Mg':
    R_per_target=CS_coherent(T,Enu,12,12)*Flux_Muller_Huber_Vogel(Enu)  #s⁻¹MeV⁻²
    R=(6.02*10**23)*R_per_target#s⁻¹MeV⁻²mol⁻¹
    return R

  elif atom =='H':
    R_per_target=CS_coherent(T,Enu,1,1)*Flux_Muller_Huber_Vogel(Enu)  #s⁻¹MeV⁻²
    R=(6.02*10**23)*R_per_target#s⁻¹MeV⁻²mol⁻¹
    return R

  elif atom =='Na':
    R_per_target=CS_coherent(T,Enu,11,12)*Flux_Muller_Huber_Vogel(Enu)  #s⁻¹MeV⁻²
    R=(6.02*10**23)*R_per_target#s⁻¹MeV⁻²mol⁻¹
    return R
  elif atom =='Cl':
    R_per_target=CS_coherent(T,Enu,17,18)*Flux_Muller_Huber_Vogel(Enu)  #s⁻¹MeV⁻²
    R=(6.02*10**23)*R_per_target#s⁻¹MeV⁻²mol⁻¹
    return R

  else:
    print("Invalid atom!")
    return 0
     

def Recoil_spectrum_atom(T,atom):
  Recoil=integrate.quad(lambda Enu: Event_rate_atom(Enu,T,atom),0,8,epsabs=1.49e-20) #s⁻¹MeV⁻¹g⁻¹
  #print(Recoil)
  return Recoil #s⁻¹MeV⁻¹mol⁻¹


##################################### Background #############################################

#### nu-bkg (Solar, Atm, DSNB) ####
DMU.loadNeutrinoFlux()
def Event_rate_nu_background(Enu,T,mineral,source):
  fluxID = DMU.nu_source_list[source]
  R_per_target=CS_coherent_mineral(T,Enu,mineral)*DMU.neutrino_flux_list[fluxID](Enu)  #s⁻¹MeV⁻²
  if mineral =='SiO2':
    n_targets_density= (6.02*10**23)/SiO2_molar_mass #SiO2 targets/g
  elif mineral =='MgSO4_7_H2O':
    n_targets_density= (6.02*10**23)/MgSO4_7_H2O_molar_mass #MgSO4 targets/g
  elif mineral =='NaCl':
    n_targets_density= (6.02*10**23)/NaCl_molar_mass #MgSO4 targets/g
  else:
    print("Invalid Mineral!")
    return 0

  R=n_targets_density*R_per_target#s⁻¹MeV⁻²g⁻¹
  return R

def Recoil_spectrum_nu_background(T,mineral,source):
  fluxID = DMU.nu_source_list[source]
  Emin,Emax=DMU.Enu_min[fluxID], DMU.Enu_max[fluxID]
  Recoil=integrate.quad(lambda Enu: Event_rate_nu_background(Enu,T,mineral,source),Emin,Emax,epsabs=1.49e-20) #s⁻¹MeV⁻¹g⁻¹
  #print(Recoil)
  return Recoil #s⁻¹MeV⁻¹g⁻¹

def Recoil_spectrum_nu_background_total(T,mineral):
  nu_sources=["8B", "hep", "atm", "DSNB", "15O", "17F", "13N", "pp"]
  R_tot=0
  R_tot_error=0
  for source in nu_sources:
    R_tot_aux=Recoil_spectrum_nu_background(T,mineral,source)
    R_tot=R_tot_aux[0]+R_tot
    R_tot_error=R_tot_aux[1]+R_tot_error
  return R_tot, R_tot_error
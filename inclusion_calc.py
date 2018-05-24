import numpy as np
import pandas as pd
import pdb
from scipy.optimize import fsolve


filename = "thermodynamic_properties.xlsx"
sheetname = "active_sheet_abbrevations"

P_amb = 0.001
R = 8.3144598
T_amb = 298.15
aT = -1.884705*1e-3
aTH = 0.706630*1e-3
beta = 0.5
beta_high = 0.25

def read_therm_table(filename, sheet):
    df = pd.read_excel(filename, sheet)
    return df

def func_k0_T(k0, dk_dT, T):
    k0_T = k0 + dk_dT * (T - 298.15)
    return k0_T

def eos_consts(k0_T, k0_1, k0_2):
    a = (1.0 + k0_1) / (1 + k0_1 + k0_T * k0_2)
    b = k0_1 / k0_T - k0_2 / (1.0 + k0_1)
    c = (1.0 + k0_1 + k0_T * k0_2) / (k0_1** 2 + k0_1 - k0_T * k0_2)
    return (a, b, c)

def lambda_nonlin(x, n, sf, T, delta_h, w):
    obj_func = (sf*(n / (n + 1))*R*T
                * (np.log(n - n*x) + np.log(1.0 - x) - np.log(1.0 + n*x)
                   - np.log(n + x)) + w*(2*x - 1.0) + delta_h)
    return obj_func

class BarInclusion(object):

    def __init__(self):
        self.df = read_therm_table(filename, sheetname)
        self.valid_phases = self.df['abbreviations'].values
        self.phase_aliases = self.df['End-member'].apply(lambda x: x.lower()).values
        self.phase_convert = {x: y for x, y in zip(self.phase_aliases, self.valid_phases)}
        self.calc_master = {}
        self.shear_master = {}
        self.T_calc_master = {}
        self.mix_master = {}

    def phase_std(self, phase):
        phase = phase.lower()
        if phase in self.valid_phases:
            return phase
        elif phase in self.phase_aliases:
            return self.phase_convert[phase]
        else:
            return 'Error'

    def hp_teos(self, P, T, phase):
        P = np.float(P)
        T = np.float(T)
        phase_proper = self.phase_std(phase)
        if phase_proper == 'Error':
            raise KeyError("That phase does not exist in the database!")

        if phase_proper not in self.calc_master.keys():
            self.calc_master[phase_proper] = {}

        if (P, T) not in self.calc_master[phase_proper].keys():
            df_p = self.df.loc[self.df['abbreviations'] == phase_proper, :].reset_index(drop=True)
            k0_T = func_k0_T(df_p['k0'], df_p['dK_dT'], T)
            a, b, c = eos_consts(k0_T, df_p['k0_1'], -df_p['k0_1']/df_p['k0'])
            if df_p['thermal_eos'][0] == 1:
                einstein_temperature = 10636.0 / (df_p['s0'] / df_p['sum_apfu'] + 6.44)
                u = einstein_temperature / T
                u_ambient = einstein_temperature / 298.15
                einstein_function_ambient = (u_ambient**2 * np.exp(u_ambient)
                                             / (np.exp(u_ambient) - 1.0)**2)
                thermal_pressure = ((df_p['alpha0']*1e-5 * k0_T
                                     * (einstein_temperature / einstein_function_ambient))
                                    * (1.0 / (np.exp(u) - 1.0) - 1.0 / (np.exp(u_ambient) - 1.0)))
                molar_volume = (df_p['v0']
                                * (1.0 - a * (1.0 - (1.0 + b * (P - thermal_pressure))**(-c))))
                if df_p['lambda'][0] == 1:
                    crit_temp_P = df_p['crit_temp'] + (df_p['max_v'] / df_p['max_s']/1000) * P
                    Q2_ambient = np.sqrt(1 - 298.15 / df_p['crit_temp'])
                    if T > crit_temp_P[0]:
                        Q2 = 0
                    else:
                        Q2 = np.sqrt((crit_temp_P - T) / df_p['crit_temp'])
                    excess_v = df_p['max_v'] * (Q2_ambient - Q2)
                    molar_volume = molar_volume + excess_v

                elif df_p['lambda'][0] == 2:
                    obj_func = lambda x: lambda_nonlin(x, df_p['n'], df_p['sf'], T,
                                                       df_p['delta_h'], df_p['w'])
                    Q = fsolve(obj_func, 0.001)
                    disorder_v = (1 - Q) * df_p['delta_v'] + df_p['wv']*Q*(1.0 - Q)
                    molar_volume = molar_volume + disorder_v
            elif df_p['thermal_eos'][0] == 2:
                molar_volume_thermal = df_p['v0'] * (1.0 + df_p['alpha0']*1e-5 * (T - 298.15))
                molar_volume = molar_volume_thermal * (1.0 - a * (1.0 - (1.0 + b * P)** (-c)))
                if df_p['End-member'] == 'Quartz':
                    pgpa = P / 10
                    crit_temp_P = df_p['crit_temp'] + 270.0 * pgpa - 10.0 * pgpa**2
                    if T > crit_temp_P[0]:
                        Ev = aTH * np.abs((crit_temp_P - T))**beta_high
                    else:
                        Ev = aT * np.abs((crit_temp_P - T))**beta
                    molar_volume = molar_volume * (1.0 + Ev)
            self.calc_master[phase_proper][(P, T)] = molar_volume[0]
        return self.calc_master[phase_proper][(P, T)]

    def Pfoot_calc(self, P_incl, P_foot, inclusion, host):
        mv_incl_T_amb_P_incl = 0
        mv_incl_T_amb_P_foot = 0
        mv_host_amb = 0
        mv_host_T_amb_P_foot = 0
        shear_modulus = 0
        for phase, mol_frac in inclusion.items():
            mv_incl_T_amb_P_incl += self.hp_teos(P_incl, T_amb, phase) * mol_frac
            mv_incl_T_amb_P_foot += self.hp_teos(P_foot, T_amb, phase) * mol_frac
        for phase, mol_frac in host.items():
            mv_host_amb += self.hp_teos(P_amb, T_amb, phase) * mol_frac
            mv_host_T_amb_P_foot += self.hp_teos(P_foot, T_amb, phase) * mol_frac
            shear_modulus += self.host_shear_calc(phase)

        mixed_objective = ((0.75/shear_modulus)*(P_incl - P_amb)
                           - (mv_incl_T_amb_P_incl / mv_incl_T_amb_P_foot
                              - mv_host_amb / mv_host_T_amb_P_foot))
        return mixed_objective



    def host_shear_calc(self, phase):
        phase_proper = self.phase_std(phase)
        if phase_proper not in self.shear_master.keys():
            self.shear_master[phase_proper] = {}
            df_p = self.df.loc[self.df['abbreviations'] == phase_proper, :].reset_index(drop=True)
            shear_modulus = df_p['shear_modulus'][0]
            poisson_ratio = df_p['poisson_ratio'][0]
            if ~np.isnan(shear_modulus):
                self.shear_master[phase_proper] = shear_modulus
            elif ~np.isnan(poisson_ratio):
                self.shear_master[phase_proper] = (3.0 * df_p['k0'][0]*(1.0 - 2.0 * poisson_ratio)) / (2.0 * (1.0 + poisson_ratio))
            else:
                self.shear_master[phase_proper] = 0.0
        return self.shear_master[phase_proper]


    def mix_elastic(self, P_incl, P_entrap, T_entrap, inclusion, host):
        # inclusion and host should be dictionaries.
        # keys are host abbreviations, values are mole fractions.
        # e.g. host = {'abv': 0.01, 'z': 0.99}
        # e.g. inclusion = {'q': 1}
        obj_func = lambda x: self.Pfoot_calc(P_incl, x, inclusion, host)
        P_footy = fsolve(obj_func, 0.001)
        mv_incl_T_amb_P_incl = 0
        mv_incl_T_entrap_P_entrap = 0
        mv_incl_T_amb_P_footy = 0
        mv_host_amb = 0
        mv_host_T_entrap_P_entrap = 0
        shear_modulus = 0
        for phase, mol_frac in inclusion.items():
            mv_incl_T_amb_P_incl += self.hp_teos(P_incl, T_amb, phase)*mol_frac
            mv_incl_T_entrap_P_entrap += self.hp_teos(P_entrap, T_entrap,  phase)*mol_frac
            mv_incl_T_amb_P_footy += self.hp_teos(P_footy, T_amb, phase)*mol_frac
        for phase, mol_frac in host.items():
            mv_host_amb += self.hp_teos(P_amb, T_amb, phase)*mol_frac
            mv_host_T_entrap_P_entrap += self.hp_teos(P_entrap, T_entrap, phase)*mol_frac
            shear_modulus += self.host_shear_calc(phase)

        mixed_objective = ((mv_incl_T_amb_P_incl / mv_incl_T_entrap_P_entrap
                            - mv_host_amb / mv_host_T_entrap_P_entrap)
                           * (mv_incl_T_entrap_P_entrap / mv_incl_T_amb_P_footy)
                            - (0.75 / shear_modulus) * (P_incl - P_amb))
        return mixed_objective


    def P_entrap_calc(self, P_incl, T_entrap, inclusion, host):
        obj_func = lambda x: self.mix_elastic(P_incl, x, T_entrap, inclusion, host)
        P_entrapment = fsolve(obj_func, 0.001)[0]
        return P_entrapment

    def T_entrap_calc(self, P_incl, P_entrap, inclusion, host):
        obj_func = lambda x: self.mix_elastic(P_incl, P_entrap, x, inclusion, host)
        T_entrapment = fsolve(obj_func, 273)[0]
        return T_entrapment


    def isobar_calc(self, T_range, P_range, inclusion, host):
        isobar_results = []
        T_use, P_use = np.meshgrid(T_range, P_range)
        for T_entrap, P_incl in zip(T_use.flatten(), P_use.flatten()):
            isobar_results.append(self.P_entrap_calc(P_incl, T_entrap, inclusion, host))
        return np.asarray(isobar_results).reshape(np.shape(T_use)), T_use, P_use









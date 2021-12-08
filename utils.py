import numpy as np
import pandas as pd
from n2o_nerp_tables import (n2o_fuel_use_df,
                             n2o_emissions_per_hectare_varying_increase_in_yield,
                             n2o_ghg_coeff_tranporting_handling,
                             n2o_crop_residue_factors_from_holos_methodology,
                             ecodistrict_factors_for_alberta,
                             reduction_modifier_df
                             )


class PROJECT_SAMPLE:
    def __init__(self,
                 area: float = 0.0,
                 synthetic_n: float = 0.0,
                 manure: float = 0.0,
                 crop: str = '',
                 ecodistrict: int = 0,
                 yield_per_acre: float = 0.0, 
                 reduction_modifier: float = 0.85):
        self._area = area
        self._synthetic_n = synthetic_n
        self._manure = manure
        self._crop = crop
        self._ecodistrict = ecodistrict
        self._yield_per_acre = yield_per_acre
        self._reduction_modifier = reduction_modifier
        self._ef = ecodistrict_factors_for_alberta[ecodistrict_factors_for_alberta.ECO_DISTRICT == self._ecodistrict]['EF_FOR_SOIL'].values.tolist()[0]/100
        self._frac_l = ecodistrict_factors_for_alberta[ecodistrict_factors_for_alberta.ECO_DISTRICT == self._ecodistrict]['FRACTION_OF_NITROGEN_LOST_IN_LEACHATE'].values.tolist()[0]/100
        self._crop_slice = n2o_crop_residue_factors_from_holos_methodology[n2o_crop_residue_factors_from_holos_methodology.CROP == self._crop]
        self._frac_of_total_dry_matter_harvested = self._crop_slice.N_FRAC_TOTAL_DRY_MATTER_HRV.values.tolist()[0]
        self._ratio_above_ground_residue_dry_matter = self._crop_slice.N_RATIO_ABOVE_GROUND_RES.values.tolist()[0]
        self._n_content_of_above_ground_residues = self._crop_slice.N_OF_ABOVE_GROUND_RES.values.tolist()[0]
        self._ratio_of_below_ground_residue_dry_matter = self._crop_slice.N_RATIO_BELOW_GROUND_RES.values.tolist()[0]
        self._n_conent_of_below_ground_residues = self._crop_slice.N_OF_BELOW_GROUND_RES.values.tolist()[0]

        self._yield = 0.0
        self._n_res = 0.0
        self._emission_manure = 0.0
        self._emission_synthetic = 0.0
        self._emission_manure_soil_crop_dynamics = 0.0
        self._n2o_nres = 0.0
        self._n_vd = 0.0
        self._n_l = 0.0
        self._emission_synthetic_soil_n_crop_dyn = 0.0
        self._total_n2o_emissions = 0.0
        self._co2e = 0.0
        self._emission_intensity = 0.0

    # Getters
    def get_area(self):
        return self._area
    def get_synthetic_n(self):
        return self._synthetic_n
    def get_manure(self):
        return self._manure
    def get_crop(self):
        return self._crop
    def get_ecodistrict(self):
        return self._ecodistrict
    def get_yield_per_acre(self):
        return self._yield_per_acre
    def get_reduction_modifier(self):
        return self._reduction_modifier
    def get_emission_factor(self):
        return self._ef
    def get_frac_l(self):
        return self._frac_l
    def get_calculated_yield(self):
        return self._yield
    def get_n_res(self):
        return self._n_res
    def get_emissions_from_manure(self):
        return self._emission_manure
    def get_emissions_synthetic_N_fertilizers(self):
        return self._emission_synthetic
    def get_emission_manure_soil_crop_dynamics(self): #
        return self._emission_manure_soil_crop_dynamics
    def get_n2o_nres(self):
        return self._n2o_nres
    def get_n_vd(self):
        return self._n_vd
    def get_n_l(self):
        return self._n_l
    def get_emission_synthetic_soil_n_crop_dyn(self):
        return self._emission_synthetic_soil_n_crop_dyn
    def get_total_n2o_emissions(self):
        return self._total_n2o_emissions
    def get_co2e(self):
        return self._co2e
    def get_emission_intensity(self):
        return self._emission_intensity


    # Setters
    def set_area(self, area):
        self._area = area
    def set_synthetic_n(self, synthetic_n):
        self._synthetic_n = synthetic_n
    def set_manure(self, manure):
        self._manure = manure
    def set_crop(self, crop):
        self._crop = crop
    def set_ecodistrict(self, ecodistrict):
        self._ecodistrict = ecodistrict
    def set_yield_per_acre(self, yield_per_acre):
        self._yield_per_acre = yield_per_acre
    def set_reduction_modifier(self, reduction_modifier):
        self._reduction_modifier = reduction_modifier
    def set_emission_factor(self, ef):
        self._ef = ef
    def set_frac_l(self, leach_frac):
        self._frac_l = leach_frac
        
    
    # Calculations
    def calculate_yield(self):
        self._yield = self._yield_per_acre * 50 /2.2 * 2.47 * self._area
    def calculate_n_res(self):
        "This equation needs to be checked again"
        A = self._yield
        B = (1 / self._frac_of_total_dry_matter_harvested)
        C = (self._ratio_above_ground_residue_dry_matter * self._n_content_of_above_ground_residues)
        D = (self._ratio_of_below_ground_residue_dry_matter * self._n_conent_of_below_ground_residues)
        E = 1/(C + D)
        F = A * B
        self._n_res = F/E        
    def calculate_emissions_from_manures(self):
        print("Manure: ", self._manure)
        print("Emission factor: ", self._ef)
        self._emission_manure = self._manure * self._ef * (44/28)
    def calculate_emission_from_synthetic_n(self):
        self._emission_synthetic = self._synthetic_n * self._ef *(44/28)
    def calculate_emissions_manure_from_soil_crop_dynamics(self):
        self._emission_manure_soil_crop_dynamics = ((self._manure * 0.2 * 0.01) + (self._manure * self._frac_l * 0.0025))*(44/28)
    def calculate_n2o_residues(self):
        self._n2o_nres = self._n_res * self._ef
    def calculate_nitrogen_vd(self):
        self._n_vd = self._synthetic_n * 0.1
    def calculate_n_leached(self):
        self._n_l = (self._synthetic_n + self._n_res) * 0.19
    def calculate_emission_synthetic_soil_crop_dynamics(self):
        self._emission_synthetic_soil_n_crop_dyn = ((self._n2o_nres)+(self._n_vd*0.01)+(self._n_l*0.025))*(44/28)
    def calculate_total_n2o_emissions(self):
        self._total_n2o_emissions = self._emission_manure + self._emission_synthetic + self._emission_manure_soil_crop_dynamics + self._emission_synthetic_soil_n_crop_dyn
    def calculate_co2_equivalent(self):
        self._co2e = self._total_n2o_emissions * 310 /1000
    def calculate_emission_intensity(self):
        self.emission_intensity = (self._co2e*self._reduction_modifier) / self._yield
    def calculate_co2e_stocks_from_n(self):
        self.calculate_yield()
        self.calculate_n_res()
        self.calculate_emissions_from_manures()
        self.calculate_emission_from_synthetic_n()
        self.calculate_emissions_manure_from_soil_crop_dynamics()
        self.calculate_n2o_residues()
        self.calculate_nitrogen_vd()
        self.calculate_n_leached()
        self.calculate_emission_synthetic_soil_crop_dynamics()
        self.calculate_total_n2o_emissions()
        self.calculate_co2_equivalent()
        self.calculate_emission_intensity()

if __name__ == "__main__":
    ps = PROJECT_SAMPLE(area=0.4,
                        synthetic_n=33.3,
                        manure=0,
                        crop='canola',
                        ecodistrict=737,
                        yield_per_acre=40)
    ps.calculate_yield()
    n2o_yield = ps.get_calculated_yield()
    print("Yield: ", n2o_yield, "\n\n")

    ps.calculate_n_res()
    n_res = ps.get_n_res()
    print("N Res: ", n_res)

    ps.calculate_emissions_from_manures()
    emission_from_manure = ps.get_emissions_from_manure()
    print("emission_from_manure: ", emission_from_manure)

    ps.calculate_emission_from_synthetic_n()
    emissions_synthetic_N_fertilizers = ps.get_emissions_synthetic_N_fertilizers()
    print("emissions_synthetic_N_fertilizers: ", emissions_synthetic_N_fertilizers)

    ps.calculate_emissions_manure_from_soil_crop_dynamics()
    emissions_manure_from_soil_crop_dynamics = ps.get_emission_manure_soil_crop_dynamics()
    print("emissions_manure_from_soil_crop_dynamics: ", emissions_manure_from_soil_crop_dynamics)

    ps.calculate_n2o_residues()
    n2o_residues = ps.get_n2o_nres()
    print("n2o_residues: ", n2o_residues)

    ps.calculate_nitrogen_vd()
    nitrogen_vd = ps.get_n_vd()
    print("nitrogen_vd: ", nitrogen_vd)

    ps.calculate_n_leached()
    n_leached = ps.get_n_l()
    print("n_leached: ", n_leached)

    ps.calculate_emission_synthetic_soil_crop_dynamics()
    emission_synthetic_soil_crop_dynamics = ps.get_emission_synthetic_soil_n_crop_dyn()
    print("emission_synthetic_soil_crop_dynamics: ", emission_synthetic_soil_crop_dynamics)

    ps.calculate_total_n2o_emissions()
    total_n2o_emissions = ps.get_total_n2o_emissions()
    print("total_n2o_emissions: ", total_n2o_emissions)

    ps.calculate_co2_equivalent()
    co2_equivalent = ps.get_co2e()
    print("co2_equivalent: ", co2_equivalent)

    ps.calculate_emission_intensity()
    emission_intensity = ps.get_emission_intensity()
    print("emission_intensity: ", emission_intensity)


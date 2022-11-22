# Written by Rebecca L. Davies
# Contact: davies@rebeccalouise.net
# Last updated 03 August 2022
# This tool is released in conjunction with the E-XQR-30 metal absorber catalog published by Davies et al. (2022a)

import json, itertools, os, csv
from typing import Any, List, Dict, Optional, Union
from astropy.cosmology import Planck18
import astropy.constants as const
import numpy as np

# The release_data folder should be one level above the script
release_data_dir = os.path.dirname(os.path.realpath(__file__)) + '/../'

# read in the rest-frame wavelengths of common transitions
ref_transition_list: List[str] = []
line_wlength_nm_list: List[float] = []
with open(release_data_dir + 'AbsorptionPathTool/transition_wavelengths.dat') as wlength_ref:
    for row in wlength_ref:
        if row[0] == '#':
            continue
        entries: List[str] = row.split()
        ref_transition_list.append(entries[0].strip())
        line_wlength_nm_list.append(float(entries[1]))

# read in a dictionary of QSOs and redshifts
qso_redshift_dict: Dict[str, float] = {}
with open(release_data_dir + 'Tables/redshifts_and_spectral_resolutions.csv','r') as qso_redshift_list:
    csv_reader = csv.reader(qso_redshift_list)
    for row in csv_reader:
        if row[0].strip() != 'Quasar':
            qso_redshift_dict[row[0]] = float(row[1])

# get the rest-frame wavelength of a particular transition
def get_transition_rest_wl(transition: str):
    return line_wlength_nm_list[ref_transition_list.index(transition)]

# calculate the relativistic velocity difference between an absorber at z = absorber_redshift and the QSO at z = qso_redshift
def calc_absorber_vel(absorber_redshift: Union[float, np.ndarray], qso_redshift: float):
    return np.log( (1 + absorber_redshift) / (1 + qso_redshift) ) * const.c.to('km/s').value

# calculate the redshift of an absorber at a velocity separation 'vel' from the QSO at z = qso_redshift
def convert_vel_to_redshift(vel: float, qso_redshift: float):
    z_arr: np.ndarray = np.arange(qso_redshift - 1, qso_redshift + 0.2, 0.001)
    vel_arr: np.ndarray = calc_absorber_vel(z_arr, qso_redshift)
    vel_arr[vel_arr > vel] = np.nan

    return z_arr[np.nanargmax(vel_arr)]

# Parent class which enables all the path length calculations
class AbsorptionPath():
    def __init__(self, config):
        
        # The adopted value of Omega matter. The routines used in this script assume a flat cosmology.
        if 'Om0' in config:
            self.Om0 = config['Om0']
        else:
            self.Om0 = Planck18.Om0

        # Bool indicating whether or not to include quasar proximity zones in the search path. 
        if 'include_proximity_zones' in config:
            self.include_proximity_zones = config['include_proximity_zones']
        else:
            self.include_proximity_zones = False

        # Float indicating the velocity threshold separating proximate and intervening absorbers (in km/s).
        if 'proximate_velocity_threshold' in config:
            self.proximate_velocity_threshold = config['proximate_velocity_threshold']
        else:
            self.proximate_velocity_threshold = 1e4

        # Bool indicating whether or not to include BAL regions in the search path. 
        if 'include_bal_regions' in config:
            self.include_bal_regions = config['include_bal_regions']
        else:
            self.include_bal_regions = False
        
        # Bool indicating whether or not to include masked regions in the search path. 
        if 'include_masked_regions' in config:
            self.include_masked_regions = config['include_masked_regions']
        else:
            self.include_masked_regions = False
        
        # List of QSO IDs to consider.
        if 'qso_idlist' in config:
            self.qso_idlist = config['qso_idlist']
        else:
            self.qso_idlist = list(qso_redshift_dict.keys())

        # list of combinations of transitions to consider when computing the absorption path. 
        self.line_combo_list: List[List[str]] = []
        
        # grid of redshift values to use when computing the absorption path length
        self.z_ref_array: np.ndarray = np.arange(1,7,1e-4)
        # Dictionary where the key is the qso_id and the value is an ndarray with the same size as 
        # self.z_ref_array, and each channel is set to either zero (absorption path DOES NOT cover 
        # this redshift) or one (absorption path DOES cover this redshift).
        self.z_covered_bool: Dict[str, np.ndarray] = {}

    def convert_z_to_dX(self, z: float):
        # Calculate the comoving distance to an object at a redshift z using Equation 1 from Davies et al. (2022). 
        # Assumes a flat cosmology.
        Lambda0: float = 1 - self.Om0
        return (2/(3 * self.Om0)) * (self.Om0 * (1 + z)**3 + Lambda0)**(1/2)

    def get_z_masked_regions(self, qso_id: str):
        # Dictionary of wavelength regions that contain BAL absorption
        bal_intervals: Dict[str, Dict[str, List[float]]] = json.load(open(release_data_dir + 'Tables/bal_wavelength_regions.json','r'))
        # Dictionary of wavelength regions that were masked due to strong skyline/telluric contamination
        regions_to_mask_dict: Dict[str, List[List[float]]] = json.load(open(release_data_dir + 'Tables/masked_wavelength_regions.json','r'))

        # Redshift of the QSO
        qso_redshift: float = qso_redshift_dict[qso_id]

        # self.z_covered_bool is a dictionary where each key is a QSO ID and each value is a numpy array. For each index i, self.z_covered_bool[qso_id][i] will either set to 1 if the redshift given by self.z_ref_array[i] is covered by the absorption search path for the relevant quasar and ion/transitions, or 0 otherwise.
        if qso_id not in self.z_covered_bool:
            self.z_covered_bool[qso_id] = np.zeros(self.z_ref_array.size)

        # Loop through self.line_combo_list, which has type List[List[transition: str]].
        # When get_z_masked_regions is called as a method on AbsorptionPathTransitions, self.line_combo_list contains a single list which is the provided list of transitions, and only one loop will be executed. In this case, the absorption path will not include wavelength regions for which ANY of the transitions lie in one of the excluded regions (e.g. BAL or masked regions).
        # When get_z_masked_regions is called as a method on AbsorptionPathPrimaryIon;
        # - if the ion is one of the doublets (MgII, FeII, CIV or SiIV), self.line_combo_list contains a single list of the two doublet lines, and only one loop will be executed
        # - if the ion is either FeII or CII, self.line_combo_list contains multiple lists representing all two-transition combinations of the transitions listed for the relevant ion in search_transitions_dict. In this case, the absorption path will only exclude bad wavelength regions for which at least one transition of ALL line combinbations lies in the relevant region.
        for line_combo in self.line_combo_list:
            # extract the rest-frame wavelengths of the transitions
            line_wavelengths: List[float] = [get_transition_rest_wl(line) for line in line_combo]
            # sort transitions in ascending order of wavelength
            line_wavelengths.sort()
            series_wlength_min: float = line_wavelengths[0]
            series_wlength_max: float = line_wavelengths[-1]
            
            # the minimum redshift where this combination of transitions detected is the redshift where the 
            # lowest wavelength line falls at the wavelength of the quasar Lya line.
            z_min: float = (1 + qso_redshift) * get_transition_rest_wl('Ly_a') / series_wlength_min - 1
            
            # The maximum redshift is:
            if self.include_proximity_zones:
                # 5000 km/s redward of the quasar redshift if proximity zones are included
                threshold_vel: float = 5000
            else:
                # 10,000 km/s blueward of the quasar redshift if proximity zones are NOT included
                threshold_vel: float = -self.proximate_velocity_threshold

            # convert this velocity separation to a redshift
            z_max: float = convert_vel_to_redshift(threshold_vel, qso_redshift)

            # this can happen for lines that fall very close to Lya, e.g. NV
            if z_min > z_max:
                continue
            
            # create a list of wavelength ranges to exclude
            exclude_wl_ranges: List[List[float]] = []
            # If not including BAL regions, add these intervals to the list of excluded regions.
            if (not self.include_bal_regions) and (qso_id in bal_intervals):
                exclude_wl_ranges.extend(list(bal_intervals[qso_id].values()))
            # If not including masked regions, exclude anything falling within these intervals.
            if (not self.include_masked_regions) and (qso_id in regions_to_mask_dict):
                exclude_wl_ranges.extend(regions_to_mask_dict[qso_id])
            
            # convert the excluded wavelength ranges into excluded redshift intervals.
            exclude_redshift_intervals: List[List[float]] = []
            for exclude_wl_range in exclude_wl_ranges:
                # for each excluded wavelength range, the minimum redshift to exclude is given by the redshift where the shortest wavelength line enters the lower end of the masked region, and the maximum redshift to exclude is given by the redshift where the longest wavelength line exits the upper end of the masked region.
                exclude_redshift_intervals.append([exclude_wl_range[0]/series_wlength_min - 1, exclude_wl_range[1]/series_wlength_max - 1])
            
            if len(exclude_redshift_intervals) > 0:
                # combine the excluded redshift intervals into a single mask
                sort_by: List = [exclude_redshift_interval[0] for exclude_redshift_interval in exclude_redshift_intervals]
                exclude_redshift_intervals = [exclude_redshift_intervals[ind] for ind in np.argsort(sort_by)]

                redshift_range: np.ndarray = np.arange(z_min,z_max,0.001)
                mask: np.ndarray = np.ones(redshift_range.size)
                covered_z_windows: List = []
                for exclude_interval in exclude_redshift_intervals:
                    if exclude_interval[0] > z_max or exclude_interval[1] < z_min:
                        continue
                    mask[np.logical_and(redshift_range >= exclude_interval[0], redshift_range < exclude_interval[1])] = 0

                while np.any(mask == 1):
                    # go through the mask, identify contiguous regions where the mask is set to one, and append them to a list of covered redshift windows.
                    masked_redshift_range: np.ndarray = redshift_range[mask==1]
                    redshift_diff: np.ndarray = masked_redshift_range[1:] - masked_redshift_range[:-1]
                    ind_zstart: int = np.nanargmin(masked_redshift_range)
                    z_start: float = masked_redshift_range[ind_zstart]
                    if ind_zstart == masked_redshift_range.size - 1 or redshift_diff.size == 1:
                        break
                    for z_end, zrange in zip(masked_redshift_range[ind_zstart+1:], redshift_diff[ind_zstart+1:]):
                        if zrange > 0.0015:
                            break
                    covered_z_windows.append([z_start,z_end])
                    # NOTE: very important to use >= and <= here, otherwise the masking doesn't work at all!
                    mask[np.logical_and(redshift_range >= z_start, redshift_range <= z_end)] = 0
                
                if len(covered_z_windows) == 0:
                    # if no windows have been added, it means that none of the region between [z_min, z_max] has been excluded, so include it all.
                    covered_z_windows = [[z_min,z_max]]

                for z_window in covered_z_windows:
                    # for each covered window, set the entries in the relevant redshift range to 1.
                    window_start_ind: int = np.nanargmin(np.abs(self.z_ref_array - z_window[0]))
                    window_end_ind: int = np.nanargmin(np.abs(self.z_ref_array - z_window[1]))
                    self.z_covered_bool[qso_id][window_start_ind : window_end_ind] = 1

            else:
                # if no regions are excluded, then set the entries across the whole queried redshift range to 1.
                window_start_ind: int = np.nanargmin(np.abs(self.z_ref_array - z_min))
                window_end_ind: int = np.nanargmin(np.abs(self.z_ref_array - z_max))
                self.z_covered_bool[qso_id][window_start_ind : window_end_ind] = 1

    # Use the self.z_covered_bool array to determine the redshift search intervals for the given quasar and ion/transitions.
    def get_z_search_intervals(self, qso_id: str, zlower: Optional[float] = None, zupper: Optional[float] = None):
        
        # Don't re-calculate self.z_covered_bool if it has already been computed.
        if qso_id not in self.z_covered_bool:
            self.get_z_masked_regions(qso_id)
        
        # 
        z_covered_bool_temp: np.ndarray = self.z_covered_bool[qso_id].copy()
        if zlower is not None:
            z_covered_bool_temp[self.z_ref_array < zlower] = 0
        if zupper is not None:
            z_covered_bool_temp[self.z_ref_array > zupper] = 0
        
        # consolidate into a single set of intervals
        z_intervals: List = []
        z_ind: int = 0
        while True:
            if z_covered_bool_temp[z_ind] == 1:
                z_start: float = self.z_ref_array[z_ind]
                while z_covered_bool_temp[z_ind] == 1:
                    z_ind += 1
                    if z_ind == z_covered_bool_temp.size:
                        z_intervals.append([z_start, self.z_ref_array[z_ind - 1]])
                        break
                z_intervals.append([z_start, self.z_ref_array[z_ind - 1]])
            
            else:
                z_ind += 1
            
            if z_ind == z_covered_bool_temp.size:
                break

        return z_intervals

    # calculate the absorption path in dX covered by a particular quasar
    def calc_quasar_dX(self, qso_id: str, zrange: Optional[List[float]] = None):
        z_intervals: List[List[float]] = self.get_z_search_intervals(qso_id)

        total_absorption_path: float = 0
        for z_interval in z_intervals:
            if zrange is not None:
                z_interval_clamped: List = [max(zrange[0], z_interval[0]), min(zrange[1], z_interval[1])]
            else:
                z_interval_clamped: List = z_interval

            if z_interval_clamped[1] < z_interval_clamped[0]:
                continue

            absorption_path: float = (self.convert_z_to_dX(z_interval_clamped[1]) - self.convert_z_to_dX(z_interval_clamped[0]))
            total_absorption_path += absorption_path

        return total_absorption_path

    # calculate the total absorption path in dX
    def calc_total_dX(self):
        total_absorption_path: float = 0
        for qso_id in self.qso_idlist:
            total_absorption_path += self.calc_quasar_dX(qso_id)

        return total_absorption_path

    # calculate the absorption path in dX over a restricted redshift range.
    def calc_dX_over_zrange(self, zrange: List[float]):
        total_absorption_path: float = 0
        for qso_id in self.qso_idlist:
            total_absorption_path += self.calc_quasar_dX(qso_id, zrange)

        return total_absorption_path

    # calculate the absorption path in dz covered by a particular quasar
    def calc_quasar_dz(self, qso_id: str, zrange: Optional[List[float]] = None):
        z_intervals: List[List[float]] = self.get_z_search_intervals(qso_id)

        total_absorption_path: float = 0
        for z_interval in z_intervals:
            if zrange is not None:
                z_interval_clamped: List = [max(zrange[0], z_interval[0]), min(zrange[1], z_interval[1])]
            else:
                z_interval_clamped: List = z_interval

            if z_interval_clamped[1] < z_interval_clamped[0]:
                continue

            total_absorption_path += (z_interval_clamped[1]-z_interval_clamped[0])

        return total_absorption_path

    # calculate the total absorption path in dz
    def calc_total_dz(self):
        total_absorption_path: float = 0
        # loop over quasars and calculate path length covered by each one.
        for qso_id in self.qso_idlist:
            total_absorption_path += self.calc_quasar_dz(qso_id)

        return total_absorption_path

    # calculate the absorption path in dz over a restricted redshift range.
    def calc_dz_over_zrange(self, zrange: List[float]):
        total_absorption_path: float = 0
        for qso_id in self.qso_idlist:
            total_absorption_path += self.calc_quasar_dz(qso_id, zrange)

    # return redshift boundaries for bins splitting the survey evenly in dX or dz.
    def split_probed_zrange(self, config = {}):
        if 'nbins' in config:
            nbins = config['nbins']
        else:
            nbins = 2
        if 'unit' in config:
            unit = config['unit']
        else:
            unit ='X'

        if 'zrange' in config:
            zlower, zupper = config['zrange']
        else:
            zlower = None
            zupper = None

        # NOTE: the sampling here needs to be SMALL, otherwise the division between the bins will not be accurate enough!
        z_tolerance = 1e-4
        dummy_z_arr = np.arange(1, 7, z_tolerance)
        absorption_path_arr = np.zeros(dummy_z_arr.size)
        for qso_id in self.qso_idlist:
            z_intervals = self.get_z_search_intervals(qso_id, zlower = zlower, zupper = zupper)
            if len(z_intervals) == 0:
                continue

            for z_interval in z_intervals:
                for ind, dummy_z in enumerate(dummy_z_arr):
                    if dummy_z > z_interval[0] - z_tolerance and dummy_z < z_interval[1]:
                        if unit in ['X', 'x']:
                            absorption_path_arr[ind] += (self.convert_z_to_dX(dummy_z + z_tolerance) - self.convert_z_to_dX(dummy_z))
                        elif unit in ['z', 'Z']:
                            absorption_path_arr[ind] += z_tolerance

        probed_zs = dummy_z_arr[np.where(absorption_path_arr > 0)]
        bin_boundaries = [probed_zs[0]]

        absorption_path_arr_normed = absorption_path_arr / np.nansum(absorption_path_arr)
        cdf = np.nancumsum(absorption_path_arr_normed)
        bin_boundaries.extend([np.interp(percentile, cdf, dummy_z_arr) for percentile in np.arange(1,nbins)/nbins])
        bin_boundaries.append(probed_zs[-1])

        path_length_weighted_mean_redshifts = [dummy_z_arr[np.nanargmin(np.abs(cdf - percentile))] for percentile in np.arange(0.5,nbins,1)/nbins]

        bin_array = [[z1, z2] for (z1, z2) in zip(bin_boundaries[:-1], bin_boundaries[1:])]

        return bin_array, path_length_weighted_mean_redshifts

    # calculate the path-length-weighted mean redshifts of bins determined by the bin_boundaries list.
    def get_path_length_weighted_mean_redshift(self, bin_boundaries, config = {}):
        if 'unit' in config:
            unit = config['unit']
        else:
            unit ='X'
        
        path_length_weighted_mean_redshifts = []
        for (zlower, zupper) in bin_boundaries:
            # NOTE: the sampling here needs to be SMALL, otherwise the division between the bins will not be accurate enough!
            z_tolerance = 1e-4
            dummy_z_arr = np.arange(1, 7, z_tolerance)
            absorption_path_arr = np.zeros(dummy_z_arr.size)
            for qso_id in self.qso_idlist:
                z_intervals = self.get_z_search_intervals(qso_id, zlower = zlower, zupper = zupper)
                if len(z_intervals) == 0:
                    continue

                for z_interval in z_intervals:
                    for ind, dummy_z in enumerate(dummy_z_arr):
                        if dummy_z > z_interval[0] - z_tolerance and dummy_z < z_interval[1]:
                            if unit in ['X', 'x']:
                                absorption_path_arr[ind] += (self.convert_z_to_dX(dummy_z + z_tolerance) - self.convert_z_to_dX(dummy_z))
                            elif unit in ['z', 'Z']:
                                absorption_path_arr[ind] += z_tolerance

            absorption_path_arr_normed = absorption_path_arr / np.nansum(absorption_path_arr)
            cdf = np.nancumsum(absorption_path_arr_normed)

            path_length_weighted_mean_redshifts.append(dummy_z_arr[np.nanargmin(np.abs(cdf - 0.5))])

        return path_length_weighted_mean_redshifts

class AbsorptionPathTransitions(AbsorptionPath):
    def __init__(self, 
            transitions: List[str], # A list of transitions to compute the absorption search path for.
            config: Dict # Configuration dictionary containing customizable parameters
        ):
        
        # initialize the parent class first
        AbsorptionPath.__init__(self, config)

        for transition in transitions:
            if transition not in ref_transition_list:
                print(f'Error! {transition} was not found in transition_wavelengths.dat. Please check that the transition is entered correctly and add an entry to the .dat file if required.')
                exit()

        self.line_combo_list = [transitions]

class AbsorptionPathPrimaryIon(AbsorptionPath):
    def __init__(self, 
            ion: str, # The ion to compute the absorption search path for.
            config: Dict # Configuration dictionary containing customizable parameters
        ):
        
        # initialize the parent class first
        AbsorptionPath.__init__(self, config)

        # transitions that were used when peforming search for absorption systems in Davies+22
        search_transitions_dict: Dict[str, str] = {
            'MgII': "MgII_2796,MgII_2803",
            'FeII': "FeII_2344,FeII_2382,FeII_2586,FeII_2600",
            'CII': "SiII_1260,OI_1302,CII_1334,SiII_1526,AlII_1670",
            'CIV': "CIV_1548,CIV_1550",
            'SiIV': "SiIV_1393,SiIV_1402",
            'NV': "NV_1238,NV_1242",
        }

        if ion not in search_transitions_dict:
            print(f'Error! {ion} is not a primary ion. Please choose one of MgII, FeII, CII, CIV, SiIV, or NV. For other ions, use AbsorptionPathTransition with a list of transitions.')
            exit()

        transition_list: List = search_transitions_dict[ion].split(',')
        if ion == 'CII':
            self.line_combo_list: List = [['CII_1334', line] for line in transition_list if line != 'CII_1334']
        elif ion in search_transitions_dict:
            self.line_combo_list: List = list(itertools.combinations(transition_list, 2))

# calls the appropriate routine depending on the desired absorption path unit and redshift range to consider.
def compute_absorption_path(absorption_path_calculator, config = {}):
    if 'unit' not in config:
        unit = 'X'

    else:
        unit = config['unit']

        if unit not in ['X','x','Z','z']:
            print(f"{unit} is not a valid unit; please select 'X' or 'z'")
            return None
    
    if 'zrange' not in config:
        if unit in ['X','x']: 
            return absorption_path_calculator.calc_total_dX()
        elif unit in ['Z','z']:
            return absorption_path_calculator.calc_total_dz()
    else:
        zrange = config['zrange']
        if unit in ['X','x']: 
            return absorption_path_calculator.calc_dX_over_zrange(zrange)
        elif unit in ['Z','z']:
            return absorption_path_calculator.calc_dz_over_zrange(zrange)

# function to compute the absorption path covered by a specific subsets of quasars from the E-XQR-30 survey for a given list of transitions.
def compute_absorption_path_for_transitions(transitions, config = {}):
    absorption_path_calculator = AbsorptionPathTransitions(transitions, config = config)

    return compute_absorption_path(absorption_path_calculator, config = config)

# function to compute the absorption path covered by a specific subsets of quasars from the E-XQR-30 survey for one of the primary ions on which the absorber catalog is based (MgII, FeII, CII, CIV, SiIV and NV).
def compute_absorption_path_for_primary_ion(ion, config = {}):
    absorption_path_calculator = AbsorptionPathPrimaryIon(ion, config = config)

    return compute_absorption_path(absorption_path_calculator, config = config)

# Functions to compute the redshift boundaries for nbins bins that split the absorption search path evenly in either dX or dz.
def split_zrange_for_transitions(transitions, config = {}):
    absorption_path_calculator = AbsorptionPathTransitions(transitions, config = config)
    
    return absorption_path_calculator.split_probed_zrange(config = config)

def split_zrange_for_primary_ion(ion, config = {}):
    absorption_path_calculator = AbsorptionPathPrimaryIon(ion, config = config)
    
    return absorption_path_calculator.split_probed_zrange(config = config)

# Functions to compute the path-length-weighted mean redshift for a set of redshift bins defined by the bin_boundaries list.
def get_mean_redshift_for_transitions(transitions, bin_boundaries, config = {}):
    # check whether list is nested
    if not any(isinstance(i, list) for i in bin_boundaries):
        bin_boundaries = [bin_boundaries]
    
    absorption_path_calculator = AbsorptionPathTransitions(transitions, config = config)
    return absorption_path_calculator.get_path_length_weighted_mean_redshift(bin_boundaries, config = config)

def get_mean_redshift_for_ion(ion, bin_boundaries, config = {}):
    # check whether list is nested
    if not any(isinstance(i, list) for i in bin_boundaries):
        bin_boundaries = [bin_boundaries]
    
    absorption_path_calculator = AbsorptionPathPrimaryIon(ion, config = config)
    return absorption_path_calculator.get_path_length_weighted_mean_redshift(bin_boundaries, config = config)
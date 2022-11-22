# Written by Rebecca L. Davies
# Contact: davies@rebeccalouise.net
# Last updated 03 August 2022
# This tool is released in conjunction with the E-XQR-30 metal absorber catalog published by Davies et al. (2022a)

# import helper scripts 
from absorption_path_tool import *

# This tool allows you to calculate the absorption path covered over any redshift range by any subset of the 42 quasars in the E-XQR-30 sample, for either:
# 1) One of the primary ions (MgII, FeII, CII, CIV, SiIV, and NV), using 'compute_absorption_path_for_primary_ion'. This is what should be used to reproduce the absorption path length values published in Table 1 of Davies et al. (2022a).
# 2) Any list of transitions, using 'compute_absorption_path_for_transitions'.

# There are many customizable options: you can choose whether or not to include quasar proximity zones, BAL regions, and masked regions in the absorption path, adjust the adopted value of Om0 (assuming a flat cosmology), and set whether the absorption path should be calculated in units of X (comoving) or z.
# To customize the behaviour, uncomment the lines of code that set the relevant optional keyword argument(s) in the 'config' dictionary below. Parameters that are not specified manually will be set to their default values indicated in the comments.
config = {
    # 'qso_idlist': ['SDSSJ0836+0054'], 
    # A list of quasar IDs to include when calculating the absorption path. 
    # If this is not set, all 42 quasars in the E-XQR-30 sample will be used.
    
    # 'zrange': [4.6, 6.2], 
    # A list of two floats indicating the lower and upper limits of the redshift range over which to calculate the absorption path. 
    # If this is not set, the absorption path will be calculated for the whole redshift range over which the relevant ion is accessible.

    # 'include_proximity_zones': False, 
    # A boolean indicating whether or not to include quasar proximity zones (with velocity offsets from the quasar redshift between -10,000 km/s and +5,000 km/s) in the absorption path.
    # The default value is False.

    'proximate_velocity_threshold': 5000.,
    # A float indicating the velocity threshold separating proximate and intervening absorbers (in km/s).
    # The default value is 10,000.

    # 'include_bal_regions': False,
    # A boolean indicating whether or not to include quasar broad absorption line (BAL) regions.
    # The default value is False.
    
    # 'include_masked_regions': False, 
    # A boolean indicating whether or not to include wavelength regions that were masked in the line search (e.g. due to strong skyline/telluric contamination).
    # The default value is False.
    
    # from astropy.cosmology import Planck18
    # 'Om0': Planck18.Om0, 
    # The value of Om0 to adopt. 
    # The default value is the Planck+2018 value (0.30966)
    # Note that the routine assumes a flat cosmology.
    
    # 'unit': 'X',
    # The unit in which to calculate the absorption path length. Can be either 'X' (default) or 'z'.

    # 'nbins': 2,
    # The number of bins over which to split the absorption path. (This parameter is only relevant when calling split_zrange.)
    # The default value is 2.
}

# Example usage for a primary ion:
absorption_path = compute_absorption_path_for_primary_ion(
        'CIV', # The ion to consider. Can be one of 'MgII','FeII','CII','CIV','SiIV', or 'NV'. 
        # The returned absorption path will represent the path length over which at least two of the considered transitions (and all required transitions) for this ion are accessible and do not fall in any of the excluded regions.

        # Configuration dictionary containing customizable parameters that can be set above.
        config = config
)
print(absorption_path)

# Example usage for a list of transitions:
absorption_path = compute_absorption_path_for_transitions(
        ['SiII_1260', 'SiII_1526'], # The list of transitions to consider. 
        # The returned absorption path will represent the length over which all of the listed transitions are accessible and do not fall in any of the excluded regions. 

        # Configuration dictionary containing customizable parameters that can be set above.
        config = config
)

# This tool also includes functions to calculate the redshift boundaries that split the absorption path evenly in either dX or dz for an arbitrary number of bins, using either split_zrange_for_transitions or split_zrange_for_primary_ion. The function returns a nested list of bin boundaries and a list of the path-length-weighted mean redshifts of those bins.
bin_boundaries, path_length_weighted_mean_redshifts = split_zrange_for_primary_ion(
        'SiIV', # The ion to consider.

        # Configuration dictionary containing customizable parameters that can be set above.
        config = config
)

# Finally, this tool includes functions to calculate the path-length-weighted mean redshift for redshift bins defined by the nested list bin_boundaries, using either get_mean_redshift_for_transitions or get_mean_redshift_for_ion.
# The function returns a list containing the mean redshifts of the bins.
path_length_weighted_mean_redshifts = get_mean_redshift_for_ion(
        'MgII', # The ion to consider.

        [[1.5,2.4], [3.4,5.7]], # a list of pairs of floats giving the redshift boundaries of the bins.
        
        # Configuration dictionary containing customizable parameters that can be set above.
        config = config
)
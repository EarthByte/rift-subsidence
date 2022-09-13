import warnings, os, io, sys
import argparse
from pathlib import Path
import contextlib
import numpy as np
import matplotlib.pyplot as plt

import pybacktrack
import RiftSubsidence as rs
import backstrip_utils as bs

def main(well_name, sea_level, dynamic_topography, cfg_file, out_dir ):
 
    Path(out_dir).mkdir(parents=True, exist_ok=True)#make output dir

    backstrip_well_filename = f'wells/{well_name}.txt'

    warnings.filterwarnings("ignore", "Well thickness .* is larger than the total sediment thickness")
    warnings.filterwarnings("ignore", ".* does not cover, or cannot interpolate, well location .*")
    
    well, decompacted_wells = pybacktrack.backstrip_well(
        backstrip_well_filename,
        pybacktrack.BUNDLE_LITHOLOGY_FILENAMES,
        sea_level_model=sea_level)

    dynamic_topography_model = pybacktrack.DynamicTopography.create_from_bundled_model(
        dynamic_topography, well.longitude, well.latitude)

    #load cfg file
    rs.read_cfg_file(cfg_file)
    rs.PLABEL = well_name #change the label name 
 
    rs.run() #run RiftSubsidence

    # Read the output from running RiftSubsidence.
    rift_subsidence_output_run = bs.read_rift_subsidence_output(well_name)

    fig = plt.figure(figsize=(8,8))
    bs.plot_tectonic_subsidence(
        fig.gca(),
        decompacted_wells,
        plot_sea_level=True,
        dynamic_topography_model=dynamic_topography_model,
        plot_subsidence_minus_dynamic_topography=True,
        rift_subsidence_output=rift_subsidence_output_run,
        offset=0)
    plt.title(f'{well_name}')
    fig.savefig(f'{out_dir}/{well_name}_residual.png')
    print(f'The figure has been saved to {out_dir}/{well_name}_residual.png')
    
if __name__ == "__main__":
    __description__ = \
    """
        For example: python plot_tectonic_subsidence.py sunrise Haq87_SealevelCurve_Longterm M6 parameters.json my_results
    """
        
    parser = argparse.ArgumentParser(description = __description__)
    
    parser.add_argument("well_name", type=str,
                    help="well name")
    parser.add_argument("sea_level_model", type=str,
                    help="sea level model name, such as Haq87_SealevelCurve_Longterm")
    parser.add_argument("dynamic_topography_model", type=str,
                    help="dynamic topography model name, such as M6")
    parser.add_argument("cfg_file", type=str,
                    help="configuration file path")
    parser.add_argument("output_folder", type=str,
                    help="output folder name")
    args = parser.parse_args()
    
    
    main(args.well_name, args.sea_level_model, args.dynamic_topography_model, args.cfg_file, args.output_folder)
    
    
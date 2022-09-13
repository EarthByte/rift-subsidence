import warnings, os, io, sys
import argparse
from pathlib import Path
import contextlib
import numpy as np
import matplotlib.pyplot as plt

import pybacktrack
import RiftSubsidence as rs
import backstrip_utils as bs

class TectonicSubsidence:
    def __init__(self, well_name, sea_level, dynamic_topography, 
                 cfg_file, out_dir, offset=0, init_run=True):
        self.well_name = well_name
        self.sea_level = sea_level
        self.dynamic_topography = dynamic_topography
        self.cfg_file = cfg_file
        self.out_dir = out_dir
        self.offset = offset
        self.ready_to_plot = False
        
        Path(out_dir).mkdir(parents=True, exist_ok=True)#make output dir

        warnings.filterwarnings("ignore", "Well thickness .* is larger than the total sediment thickness")
        warnings.filterwarnings("ignore", ".* does not cover, or cannot interpolate, well location .*")
    
        backstrip_well_filename = f'wells/{well_name}.txt'
        self.well, self.decompacted_wells = pybacktrack.backstrip_well(
            backstrip_well_filename,
            pybacktrack.BUNDLE_LITHOLOGY_FILENAMES,
            sea_level_model=sea_level)

        self.dynamic_topography_model = pybacktrack.DynamicTopography.create_from_bundled_model(
            dynamic_topography, self.well.longitude, self.well.latitude)
        
        #load cfg file
        rs.read_cfg_file(cfg_file)
        rs.PLABEL = well_name #change the label name 
 
        if init_run:
            f = io.StringIO()
            with contextlib.redirect_stdout(f):   
                rs.run() #run RiftSubsidence

            # Read the output from running RiftSubsidence.
            self.rift_subsidence_output = bs.read_rift_subsidence_output(well_name)
            self.ready_to_plot = True

        self.BETAs = list(np.arange(rs.BETMIN,rs.BETMAX,rs.BETINC))
        
        # 'decompacted_wells' is a list of pybacktrack.DecompactedWell.
        # Extract the age and tectonic subsidence from each decompacted well in the list.
        self.ages = [decompacted_well.get_age()
                for decompacted_well in self.decompacted_wells]
        self.tectonic_subsidences = [decompacted_well.get_tectonic_subsidence()
                for decompacted_well in self.decompacted_wells]
        self.min_max_tectonic_subsidences = [decompacted_well.get_min_max_tectonic_subsidence()
                for decompacted_well in self.decompacted_wells]
        self.tectonic_subsidence_uncertainties = [0.5 * (max_tectonic_subsidence - min_tectonic_subsidence)
                for min_tectonic_subsidence, max_tectonic_subsidence in self.min_max_tectonic_subsidences]
        
        # Keep track of the dynamic topography at each age.
        #
        # NOTE: Positive values represent elevation which is the opposite of subsidence
        # (which is positive going down).
        # So we need to negate dynamic topography to turn an elevation into a depth.
        self.dynamic_topography_depths = [-self.dynamic_topography_model.sample(decompacted_well.get_age())
                for decompacted_well in self.decompacted_wells]
        
        self.rift_subsidence_model=rs
        self.well_name = well_name
        
    def run(self):
        f = io.StringIO()
        with contextlib.redirect_stdout(f):   
            rs.run() #run RiftSubsidence

        # Read the output from running RiftSubsidence.
        self.rift_subsidence_output = bs.read_rift_subsidence_output(self.well_name)
        self.ready_to_plot = True

    
    def get_rift_subsidence_model():
        return self.rift_subsidence_model
    
    
    def plot_sea_level(self,ax):
        # Keep track of the sea level (relative to present day) at each age.
        #
        # NOTE: Positive values represent sea-level rise which is the opposite of depth
        # (which is positive going down).
        # So we need to negate sea level to turn a rise into a depth.
        sea_level_depths = [-decompacted_well.get_sea_level()
                for decompacted_well in self.decompacted_wells]

        ax.plot(
            self.ages,
            sea_level_depths,
            color='blue',
            label='sea level',
            linestyle='-',
            linewidth=2.0)

        
    def plot_tectonic_subsidence(self,ax):
        ax.errorbar(
            self.ages,
            self.tectonic_subsidences,
            yerr=self.tectonic_subsidence_uncertainties,
            fmt='-o',
            color='black',
            label='tectonic subsidence',
            linestyle='-',
            linewidth=2.0)

        
    def plot_dynamic_topography(self,ax):
                  
        ax.plot(
            self.ages,
            self.dynamic_topography_depths,
            color='red',
            label='dynamic topography',
            linestyle='-',
            linewidth=2.0)
        
        
    def plot_subsidence_minus_dynamic_topography(self,ax):
        ax.errorbar(
            self.ages,
            [(self.tectonic_subsidences[age_index] - self.dynamic_topography_depths[age_index])
                 for age_index in range(len(self.ages))],
            yerr=self.tectonic_subsidence_uncertainties,
            fmt='-o',
            color='magenta',
            label='tectonic subsidence rel. dynamic topography',
            linestyle='-',
            linewidth=2.0)
        
        
    def plot_rift_subsidence(self,ax,index):
        #check if the model is ready to plot. if not, call run() first.
        if not self.ready_to_plot:
            run()
        rift_subsidence_ages, rift_subsidence_subsidences = self.rift_subsidence_output[index]
        # Plot tectonic subsidence output by RiftSubsidence program (for a specific BETA).
        ax.plot(
                rift_subsidence_ages,
                [ x + self.offset for x in rift_subsidence_subsidences],
                color='green',
                label='rift subsidence (BETA:{0})'.format(self.BETAs[index]),
                linestyle='-',
                linewidth=2.0)
        
        
    def get_residuals(self):
        #check if the model is ready to plot. if not, call run() first.
        if not self.ready_to_plot:
            run()
        try:
            return self.residuals
        except AttributeError:
            return bs.get_residuals(self.rift_subsidence_output, self.ages, 
                                     self.tectonic_subsidences, self.dynamic_topography_depths, 
                                     self.offset)
    
    
    def plot_residual(self,ax,index):
        #check if the model is ready to plot. if not, call run() first.
        if not self.ready_to_plot:
            run()
            
        self.residuals = self.get_residuals()
        
        # Plot tectonic subsidence minus dynamic topography minus rift subsidence.
        ax.errorbar(
                self.ages,
                self.residuals[index],
                yerr=self.tectonic_subsidence_uncertainties,
                fmt='-o',
                label=f'subsidence residual (BETA:{self.BETAs[index]})',
                linestyle='-',
                linewidth=2.0)
    
    
    def config_ax(self,ax):
        ax.invert_xaxis()
        ax.invert_yaxis()

        ax.set_ylabel('Depth (m)', fontsize=12)
        ax.set_xlabel('Age (Ma)', fontsize=12)
        ax.grid(linestyle='--',alpha=0.5)
        ax.legend(fontsize=10, loc='lower left')

        
def main(well_name, sea_level, dynamic_topography, cfg_file, out_dir):
 
    ts = TectonicSubsidence(well_name, sea_level, dynamic_topography, cfg_file, out_dir)
    fig, ax = plt.subplots(figsize=(8,8))
    
    ts.plot_sea_level(ax)
    ts.plot_tectonic_subsidence(ax)
    ts.plot_dynamic_topography(ax)
    ts.plot_subsidence_minus_dynamic_topography(ax)
    ts.plot_rift_subsidence(ax,0)
    ts.plot_residual(ax,0)
    ts.config_ax(ax)

    ax.set_title(f'{well_name}')
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
    
    
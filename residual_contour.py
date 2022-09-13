import warnings, os, io
import argparse
from pathlib import Path
import contextlib
import numpy as np
import matplotlib.pyplot as plt

import pybacktrack
import RiftSubsidence as rs
import backstrip_utils as bs
from tectonic_subsidence import TectonicSubsidence

class ResidualContour():
    def __init__(self, tbeg1_max, tbeg1_min, tbeg1_inc, 
                      well_name, sea_level_model, dynamic_topography_model, 
                      cfg_file, output_folder,
                      betmin=None, betmax=None, betinc=None, 
                      tend1=None, tstop=None, offset=0):
        
        self.ts = TectonicSubsidence(well_name, sea_level_model, dynamic_topography_model, 
                                      cfg_file, output_folder, offset, init_run=False)
        self.well_name = well_name
        self.tbeg1_max = tbeg1_max
        self.tbeg1_min = tbeg1_min
        self.tbeg1_inc = tbeg1_inc
        self.betmin = betmin 
        self.betmax = betmax 
        self.betinc = betinc 
        self.tend1 = tend1 
        self.tstop = tstop 
        self.offset = offset
        
        self.TBEG1s=np.arange(tbeg1_max,tbeg1_min, -tbeg1_inc)
        
        self.output_folder=output_folder
        Path(output_folder).mkdir(parents=True, exist_ok=True)#make output dir
        
        if betmin is not None: 
            self.ts.rift_subsidence_model.BETMIN = betmin 
        if betmax is not None:
            self.ts.rift_subsidence_model.BETMAX = betmax 
        if betinc is not None:
            self.ts.rift_subsidence_model.BETINC = betinc 

        if tend1 is not None:
            self.ts.rift_subsidence_model.TEND1 = tend1 
        if tstop is not None:
            self.ts.rift_subsidence_model.TSTOP = tstop 

        self.BETAs=np.arange(self.ts.rift_subsidence_model.BETMIN,
                             self.ts.rift_subsidence_model.BETMAX,
                             self.ts.rift_subsidence_model.BETINC)
        print("BETAs: ",self.BETAs)
        print("TBEG1s: ",self.TBEG1s)

        
    def run(self):
        self.data=[]
        for TBEG1 in self.TBEG1s: #for each begin time, run the RiftSubsidence program
            self.ts.rift_subsidence_model.TBEG1  = TBEG1
            f = io.StringIO()
            with contextlib.redirect_stdout(f):   
                self.ts.run() #run RiftSubsidence

            residuals = self.ts.get_residuals()
            row=[]
            for r in residuals:
                #row.append(np.mean(np.abs(r)))
                row.append(np.abs(np.mean(r)))
            self.data.append(row) 
            
            
    def plot_residual_contour(self, ax):
        #convert the result to numpy array   
        rr=np.array(self.data)
        print('The result shape: ',rr.shape)

        #plot the contour
        durations=self.TBEG1s - self.ts.rift_subsidence_model.TEND1
        print('durations:' ,durations)

        xv, yv = np.meshgrid(self.BETAs, durations)

        origin='lower'

        levels = [0, 20,50,70, 100, 150, 200, 300, 400, 500, 600, 700, 800]
        CS3 = ax.contourf(xv, yv, rr, levels,
                            cmap="gray",
                            origin=origin,
                            extend='both')

        CS4 = ax.contour(xv, yv, rr, levels,
                          colors=('r',),
                          linewidths=(1,),
                          origin=origin)

        ax.set_title('Sensitivity Analysis')
        ax.clabel(CS4, fmt='%d', colors='y', fontsize=14)
        ax.set_ylabel('Duration (TBEG1 - TEND1)')
        ax.set_xlabel('BETA')
        return CS3
    
    def save_residual_contour_plot(self):
        fig, ax = plt.subplots(figsize=(12,8))

        CS3 = self.plot_residual_contour(ax)
        fig.colorbar(CS3)
        fig.savefig(f'{self.output_folder}/{self.well_name}_residual_contour.png')
        print(f'The figure has been saved to {self.output_folder}/{self.well_name}_residual_contour.png')
    

if __name__ == "__main__":
    __description__ = \
    """
        For example: python residual_contour.py 196 90 7 sunrise Haq87_SealevelCurve_Longterm M6 parameters.json my_results
    """
        
    parser = argparse.ArgumentParser(description = __description__)
    
    parser.add_argument("tbeg1_max", type=int,
                    help="TBEG1 max")
    parser.add_argument("tbeg1_min", type=int,
                    help="TBEG1 min")
    parser.add_argument("tbeg1_inc", type=int,
                    help="TBEG1 increment")
    
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

    parser.add_argument("--betmin", type=float, help="BETMIN", default=1.)
    parser.add_argument("--betmax", type=float, help="BETMAX", default=1.301)
    parser.add_argument("--betinc", type=float, help="BETINC", default=0.03)
    parser.add_argument("--tend1", type=int, help="TEND1", default=90.)
    parser.add_argument("--tstop", type=int, help="TSTOP", default=90.)
    
    parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")
    
    args = parser.parse_args()
    
    rc = ResidualContour(args.tbeg1_max, args.tbeg1_min, args.tbeg1_inc, args.well_name, args.sea_level_model, 
         args.dynamic_topography_model, args.cfg_file, args.output_folder,
         args.betmin, args.betmax, args.betinc, args.tend1, args.tstop)
    rc.run()
    rc.save_residual_contour_plot()
    
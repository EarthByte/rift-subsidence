import numpy as np

def read_rift_subsidence_output(rift_subsidence_output_label):
    
    # Read all lines in the RiftSubsidence output file.
    rift_subsidence_output_filename = '{0}.dat'.format(rift_subsidence_output_label)
    with open(rift_subsidence_output_filename, 'r') as rift_subsidence_output_file:
        lines = rift_subsidence_output_file.readlines()
    
    beta_groups = []
    
    # First beta group starts with age and subsidence lists.
    beta_groups.append(([], []))
    
    for line in lines:
        line = line.strip()
        # If encounter empty line then we're finished.
        if not line:
            break
        
        # Get age and subsidence.
        age, subsidence, _, _ = [float(item) for item in line.split()]
        if age < 0.0:
            # Start a new beta group with empty age and subsidence lists.
            beta_groups.append(([], []))
            continue
        
        # Convert subsidence from Km to metres.
        subsidence *= 1000
        
        # Add to end of current beta group (last group in beta list).
        beta_groups[-1][0].append(age)
        beta_groups[-1][1].append(subsidence)
    
    return beta_groups

def get_residuals(rift_subsidence_output, ages, tectonic_subsidences, dynamic_topography_depths, offset):
    residuals=[]
    for beta_index, (rift_subsidence_ages, rift_subsidence_subsidences) in enumerate(rift_subsidence_output):
        interpolated_rift_data = np.interp(ages,np.flip(rift_subsidence_ages),np.flip(rift_subsidence_subsidences))
        residual = np.array(tectonic_subsidences)-np.array(dynamic_topography_depths)-np.array(interpolated_rift_data)-offset
        residuals.append(residual.tolist())
    return residuals

def get_residuals_from_wells(rift_subsidence_output, decompacted_wells, dynamic_topography_model, offset=0):
    ages = [decompacted_well.get_age()
            for decompacted_well in decompacted_wells]
    tectonic_subsidences = [decompacted_well.get_tectonic_subsidence()
            for decompacted_well in decompacted_wells]  
    dynamic_topography_depths = [-dynamic_topography_model.sample(decompacted_well.get_age())
                for decompacted_well in decompacted_wells]
    return get_residuals(rift_subsidence_output, ages, tectonic_subsidences, dynamic_topography_depths, offset)

def plot_tectonic_subsidence(
        ax,
        decompacted_wells,
        dynamic_topography_model=None,
        rift_subsidence_output=None,
        offset=0,
        plot_sea_level=False,
        plot_subsidence_minus_dynamic_topography=False,
        plot_rift_subsidence_output=True,
        plot_dynamic_topography=True,
        plot_tectonic_minus_dynamic_topography_minus_rift=False,
        plot_tectonic_subsidence_data=True):

    # 'decompacted_wells' is a list of pybacktrack.DecompactedWell.
    # Extract the age and tectonic subsidence from each decompacted well in the list.
    ages = [decompacted_well.get_age()
            for decompacted_well in decompacted_wells]
    tectonic_subsidences = [decompacted_well.get_tectonic_subsidence()
            for decompacted_well in decompacted_wells]
    min_max_tectonic_subsidences = [decompacted_well.get_min_max_tectonic_subsidence()
            for decompacted_well in decompacted_wells]
    tectonic_subsidence_uncertainties = [0.5 * (max_tectonic_subsidence - min_tectonic_subsidence)
            for min_tectonic_subsidence, max_tectonic_subsidence in min_max_tectonic_subsidences]
    
    if plot_sea_level:
        # Keep track of the sea level (relative to present day) at each age.
        #
        # NOTE: Positive values represent sea-level rise which is the opposite of depth
        # (which is positive going down).
        # So we need to negate sea level to turn a rise into a depth.
        sea_level_depths = [-decompacted_well.get_sea_level()
                for decompacted_well in decompacted_wells]
    
    if dynamic_topography_model:
        # Keep track of the dynamic topography at each age.
        #
        # NOTE: Positive values represent elevation which is the opposite of subsidence
        # (which is positive going down).
        # So we need to negate dynamic topography to turn an elevation into a depth.
        dynamic_topography_depths = [-dynamic_topography_model.sample(decompacted_well.get_age())
                for decompacted_well in decompacted_wells]


    if plot_tectonic_subsidence_data:
        # Plot tectonic subsidence.
        ax.errorbar(
            ages,
            tectonic_subsidences,
            yerr=tectonic_subsidence_uncertainties,
            fmt='-o',
            color='black',
            label='tectonic subsidence',
            linestyle='-',
            linewidth=2.0)

    if plot_sea_level:
        # Plot sea level.
        ax.plot(
            ages,
            sea_level_depths,
            color='blue',
            label='sea level',
            linestyle='-',
            linewidth=2.0)

    if dynamic_topography_model and plot_dynamic_topography:
        # Plot dynamic topography.
        ax.plot(
            ages,
            dynamic_topography_depths,
            color='red',
            label='dynamic topography',
            linestyle='-',
            linewidth=2.0)
    
    if plot_subsidence_minus_dynamic_topography:
        # Plot tectonic subsidence minus dynamic topography.
        ax.errorbar(
            ages,
            [(tectonic_subsidences[age_index] - dynamic_topography_depths[age_index])
                 for age_index in range(len(ages))],
            yerr=tectonic_subsidence_uncertainties,
            fmt='-o',
            color='magenta',
            label='tectonic subsidence rel. dynamic topography',
            linestyle='-',
            linewidth=2.0)
    
    if rift_subsidence_output and plot_rift_subsidence_output:
        for beta_index, (rift_subsidence_ages, rift_subsidence_subsidences) in enumerate(rift_subsidence_output):
            # Plot tectonic subsidence output by RiftSubsidence program (for a specific BETA).
            ax.plot(
                rift_subsidence_ages,
                [ x + offset for x in rift_subsidence_subsidences],
                color='green',
                label='rift subsidence {0}'.format(beta_index+1),
                linestyle='-',
                linewidth=2.0)
            
    if plot_tectonic_minus_dynamic_topography_minus_rift:
        colours=['red', 'blue', 'green']
        residuals=get_residuals(rift_subsidence_output,ages,tectonic_subsidences,dynamic_topography_depths,offset)
        assert len(colours) == len(residuals) #colours and residuals must have the same length
        # Plot tectonic subsidence minus dynamic topography minus rift subsidence.
        for colour, residual, beta_index in zip(colours,residuals,range(len(colours))):
            ax.errorbar(
                ages,
                residual,
                yerr=tectonic_subsidence_uncertainties,
                fmt='-o',
                color=colour,
                label=f'subsidence residual {beta_index+1}',
                linestyle='-',
                linewidth=2.0)
        ax.set_title(f'Subsidence Residual(offset:{offset})')

    ax.invert_xaxis()
    ax.invert_yaxis()

    ax.set_ylabel('Depth (m)', fontsize=12)
    ax.set_xlabel('Age (Ma)', fontsize=12)
    ax.grid(linestyle='--',alpha=0.5)
    ax.legend(fontsize=10, loc='lower left')

    return 



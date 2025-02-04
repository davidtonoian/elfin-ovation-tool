import numpy as np
import pandas as pd
import pyspedas 

#for loading subdirectory ovation data
import os
#to ignore DeprecationWarning in pyspedas.tplot from printing secondary axis
import warnings

# from astropy.constants import R_earth
# import astropy.units as u
#RE = R_earth.to(u.km).to_value() #
RE = 6378.1 #km

def load_example(i = 0, nflux_l_limit = 9):
    """
    Loads example of OVATION-prime date (from local example files) and corresponding
    ELFIN data (using pyspedas) into tplot variables (both separated and combined).
    Tplot variables can by used by executing pyspedas.pytplot.tplot([TPLOT_VARIABLE_LIST]) 
    or by using summary_plot() function
    
    Inputs:
    i: int, optional
        Number of example (from 0 to 22) 
    nflux_l_limit: int, optional
        limit of the nflux that will be turned to NaN (due to being too small)
    Returns:
    list
        List of tplot variables created.
    """
    module_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    ovation_file = f"ovation_data/OVATION_{i}.txt"
    ovation = np.loadtxt(os.path.join(module_dir, ovation_file))

    ovation_energy = ovation[0, 4:] 
    ovation = ovation[1:, :]
    ovation_time = ovation[:, 0]
    ovation_nflux = ovation[:, 4:] 
    ovation_dict = {'x': ovation_time, 
                    'y':ovation_nflux*1e3, #+switching from 1/keV to 1/MeV
                    'v':ovation_energy       
                    } 
    pyspedas.store_data('OVATION_nflux_para', data = ovation_dict,)
    pyspedas.options('OVATION_nflux_para', 
                     opt_dict ={'spec': True, 'ylog': True, 'y_range': [2.5e-2, 55],
                                'zlog': True,'z_range': [10,np.nanmax(ovation_nflux*1000)],
                                'ztitle': '#/(s-cm$^2$\n-str-MeV)',
                                'ytitle':'OVATION\nnflux_para',
                                'ysubtitle': '[keV]'}
                    )

    #download elfin data corresponding to OVATION file
    trange  = [ovation_time[0]-1, ovation_time[-1]+1] #1sec gaps added for correct download of ELFIN data
    tvars = pyspedas.projects.elfin.epd(trange=trange, probe='a', level = 'l2', fullspin = True)
    tvars_state = pyspedas.projects.elfin.state(trange=trange, probe='a')
    tvars += tvars_state
    
    #merging ELFIN&OVATION spectra
    elfin_t, elfin_y, elfin_en = pyspedas.get_data('ela_pef_fs_nflux_para')
    ovation_t, ovation_y, ovation_en = pyspedas.get_data('OVATION_nflux_para')
    #Turn small nflux into NaNs
    ovation_y = np.where(ovation_y<9, np.nan, ovation_y) 

    #store combined spectra in tplot variable
    pyspedas.store_data('ovation_ela_nflux_para', 
                        data = {'x': elfin_t, 'y': np.concatenate((ovation_y, elfin_y), axis = 1), 
                                'v': np.concatenate((ovation_en, elfin_en))}
                        )
    pyspedas.options('ovation_ela_nflux_para', 
                    opt_dict ={'spec': True, 'ylog': True, 'y_range': [2.5e-2, 6800],
                                'zlog': True,'z_range': [10,np.nanmax(ovation_nflux*1000)],
                                'ztitle': '#/(s-cm$^2$\n-str-MeV)','ytitle':'ovation_ela_\nnflux_para',
                                'ysubtitle': '[keV]'}
                    )
    tvars += ['OVATION_nflux_para', 'ovation_ela_nflux_para']

    return tvars

def summary_plot(save_png = None, xsize=6, ysize=10,**kwargs):
    """
    Creates summary plot of ELFIN+OVATION spectra assuming spectra is already
    are already loaded with load_example()
    """

    pyspedas.tplot_math.tkm2re('ela_pos_gei')
    pyspedas.split_vec('ela_pos_gei')

    #rotating from gei to sm
    pyspedas.cotrans(name_in='ela_pos_gei', name_out='ela_pos_sm',
                     coord_in='gei', coord_out='sm')
    t, y =pyspedas.get_data('ela_pos_sm')
    #Creating dipole L-shell,  MLT, MLAT (similar to elf_mlt_lat.pro)
    out = pyspedas.xyz_to_polar(y, co_latitude=True)
    r = out[:,0]
    theta = out[:, 1]
    phi = out[:, 2]
    L = r/RE/np.sin(theta/180*np.pi)**2
    mlat = 90-theta 
    MLT = 12+phi*12/180

    #storing dipole coords as tplot variables
    pyspedas.store_data('ela_pos_l_dip', data = {'x': t, 'y': L})
    pyspedas.options('ela_pos_l_dip','ytitle','L')
    pyspedas.store_data('ela_pos_mlat_dip', data = {'x': t, 'y': mlat})
    pyspedas.options('ela_pos_mlat_dip','ytitle','MLAT')
    pyspedas.store_data('ela_pos_mlt_dip', data = {'x': t, 'y': MLT})
    pyspedas.options('ela_pos_mlt_dip','ytitle','MLT')

    #to ignore DeprecationWarning in pyspedas.tplot from printing secondary axis
    #will be removed in future version
    warnings.filterwarnings('ignore', category=UserWarning)

    fig, ax =pyspedas.tplot(['ela_pef_fs_nflux_perp',
                             'ela_pef_fs_nflux_para',
                             'OVATION_nflux_para',
                             'ovation_ela_nflux_para'],
                            var_label = ['ela_pos_mlat_dip', 
                                         'ela_pos_mlt_dip', 
                                         'ela_pos_l_dip'], 
                            xsize=xsize, ysize=ysize, 
                            return_plot_objects=True,
                            **kwargs)
    ##external saving for figures (same as pyspedas but with "tight_layout" on,
    ##fixes saved images)
    if save_png:
        fig.savefig(f'{save_png}', bbox_inches = 'tight') 
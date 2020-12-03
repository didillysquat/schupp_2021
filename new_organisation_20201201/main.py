"""Starting again with the schupp project.

We now have some datasheets that should be in a more useable format.

We have a SymPortal datasheet that was used for the submission of all samples related to this project.
This is: new_organisation_20201201/20201018_schupp_all_data_datasheet.xlsx

We have the SymPortal output main directory that is: new_organisation_20201201/20201018_schupp_all

We also have the physiological data sheet that is an excel file put together by Mareen and Sam
This is: new_organisation_20201201/schupp_et_al_physiological_coral_data.xlsx

The first aim of this script will be to create two figures.
The first figure is as laid out in the online google word doc and shows the physiological data
for each of the species across the time of the study.
We will have a column for each species and a row for each of the types of data.
The species will be:
Acropora digitifera, Acropora hyacinthus, Pocillopora damicornia and Leptastrea purpurea
So 4 columns.
The rows will be:
Adult survival
Adult zooxs
Recruit survival
Recruit zooxs
Recruit growth
Recruit Fv/Fm
Not all data is available for each of the species. For the time being I think it is helpful
to leave a blank plot in these cases so that we know the data is missing.
"""

import os
import pandas as pd
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from collections import defaultdict
import numpy as np
import re
DEGREE_SIGN = u'\N{DEGREE SIGN}'


class schupp_figures:
    def __init__(self):
        # Paths
        self.root_dir = os.path.abspath(os.path.dirname(__file__))
        self.sp_datasheet_path = os.path.join(self.root_dir, '20201018_schupp_all_data_datasheet.xlsx')
        self.sp_data_path = os.path.join(self.root_dir, '20201018_schupp_all')
        self.physiological_data_path = os.path.join(self.root_dir, 'schupp_et_al_physiological_coral_data.xlsx')

        # Datatypes and species
        self.species_full = [
            'Acropora digitifera', 'Acropora hyacinthus', 'Pocillopora damicornis', 'Leptastrea purpurea'
        ]
        self.species_short = ['ad', 'ah', 'pd', 'lp']
        self.species_short_to_full_dict = {
            'ad':'Acropora digitifera', 'ah':'Acropora hyacinthus',
            'pd':'Pocillopora damicornis', 'lp':'Leptastrea purpurea'
        }
        self.data_types = ['adult_survival', 'adult_zooxs', 'recruit_survival', 'recruit_size_fv_fm', 'recruit_zooxs']

        # Dataframe
        # We want to create one single dataframe that holds all of the physiology data.
        # We will hold the zooxs data in seperate dataframes that we will create from the SP outputs
        # For the physiology data we will need a dataframe for survival and one for growth and fv_fm
        # This is because survival was assesed as a batch (i.e. 3 dead out of 5), while the growth and fv_fm
        # were done on a persample basis.

        fv_fm_size_data, survival_data = self._populate_data_holders()
        self.survival_df, self.fv_fm_size_df = self._make_survival_fv_fm_size_dfs(fv_fm_size_data, survival_data)

        # Figure
        # 10 wide, 6 deep
        self.fig = plt.figure(figsize=(10, 6))
        self.gs = gridspec.GridSpec(5, 4)
        self.axes = [[] for sp in self.species_short]
        # We will have the various axes as coordinates in columns and rows
        for i, sp in enumerate(self.species_short):
            for j, data_type in enumerate(self.data_types):
                # create a plot
                self.axes[i].append(plt.subplot(self.gs[j, i]))
        self.temp_plot_marker_dict = {29: "v", 30: "o", 31: "^", '29': "v", '30': "o", '31': "^"}
        self.temp_plot_colour_dict = {29: "#a6a6a6", 30: "#898989", 31: "#0d0d0d", '29': "#a6a6a6", '30': "#898989", '31': "#0d0d0d"}

    def plot_fig_1(self):
        """
        The main plotting method for plotting up the figure that has the four species columns with a plot
        for each of the data types.
        """
        for i, sp in enumerate(self.species_short):
            for j, data_type in enumerate(self.data_types):
                if data_type == 'adult_survival':
                    self._plot_adult_survival(i, j, sp)

    def _plot_adult_survival(self, i, j, sp):
        if sp in ['ad', 'ah']:
            # Plot the aduult_survival as just that, survival. So start with 100% and then decrease
            ax = self.axes[j][i]
            # We will need to do an ax.errorbar for each of the temperatures
            working_df = self.survival_df[
                (self.survival_df['adult_recruit'] == 'adult') & (self.survival_df['species'] == sp)]
            working_df['survival_percent'] = (working_df['survival'] / 5) * 100
            for temp in working_df['temperature'].unique():
                ser = working_df[working_df['temperature'] == temp]
                # Calc average survival for each time point and the standard error of the mean
                means = [ser[ser['time_value'] == time_val]['survival_percent'].mean() for time_val in
                         ser['time_value'].unique()]
                sem = [ser[ser['time_value'] == time_val]['survival_percent'].sem() for time_val in
                       ser['time_value'].unique()]
                ax.errorbar(
                    x=ser['time_value'].unique(), y=means, yerr=sem, marker='o',
                    linestyle='--', linewidth=1, ecolor=self.temp_plot_colour_dict[temp],
                    elinewidth=1, color=self.temp_plot_colour_dict[temp], markersize=2, label=f'{temp}{DEGREE_SIGN}C'
                )
            ax.legend(loc='lower left', fontsize='xx-small')
            ax.set_xlabel('days', fontsize='xx-small')
            if sp == 'ad':
                # only set one per row
                ax.set_ylabel('survial %', fontsize='xx-small')
            ax.set_title(self.species_short_to_full_dict[sp], fontsize='small')
            plt.tight_layout()
            foo = 'bar'
        else:
            # Then the data does not exist.
            pass

    def _make_survival_fv_fm_size_dfs(self, fv_fm_size_data, survival_data):
        survival_df = pd.DataFrame(
            data=survival_data,
            columns=['species', 'adult_recruit', 'time_value', 'time_unit', 'tank', 'temperature', 'exp_type',
                     'survival'])
        fv_fm_size_df = pd.DataFrame(
            data=fv_fm_size_data,
            columns=[
                'species', 'adult_recruit', 'time_value', 'time_unit', 'temp', 'tank', 'rack',
                'rack_row', 'rack_col', 'cyl_vol', 'fv_fm', 'exp_type',
            ])
        survival_df['survival_percent'] = 5 - survival_df['survival'] * 100
        return survival_df, fv_fm_size_df

    def _populate_data_holders(self):
        survival_data = []
        fv_fm_size_data = []
        for sp in self.species_short:
            for data_type in self.data_types:
                if data_type == 'adult_survival':
                    if sp in ['ad', 'ah']:
                        self._populate_adult_survival(sp, survival_data)
                    else:
                        # There is no adult survival data for pd or lp
                        pass

                elif data_type == 'recruit_survival':
                    self._populate_recruit_survival(sp, survival_data)

                elif data_type == 'recruit_size_fv_fm':
                    fv_fm_size_dd = self._populate_recruit_fv_fm_size_dd(sp)
                    self._populate_recruit_fv_fm_size_data(fv_fm_size_data, fv_fm_size_dd, sp)
        return fv_fm_size_data, survival_data

    def _populate_recruit_fv_fm_size_data(self, fv_fm_size_data, fv_fm_size_dd, sp):
        # At this point we have the dict populated and we can now populate the fv_fm_size_data
        for k in fv_fm_size_dd.keys():
            temp, tank, rack, rack_row, rack_col, time = k.split('_')
            fv_fm_size_data.append([
                sp, 'recruit', int(time), 'month',
                temp, tank, rack, rack_row, rack_col,
                fv_fm_size_dd[k]['cylinder_vol'], fv_fm_size_dd[k]['yield'], 'main'
            ])

    def _populate_recruit_fv_fm_size_dd(self, sp):
        # A dictionary where identifier is key and value is a dict with keys of 'fv_fm' and 'cyl_vol
        # Identifier will be made in the form "{temp}_{tank}_{rack}_{rack_row}_{rack_row}"
        # We will place numpy.np for val if no val is available
        # We will check to see that there are not duplicate entries for a given sample and data_type
        fv_fm_size_dd = defaultdict(dict)
        xl_df = pd.read_excel(io=self.physiological_data_path, sheet_name=f'{sp}_recruit_size_fv_fm_bh')
        for ind in xl_df.index:
            self._populate_fv_fm_sample_rows(ind=ind, fv_fm_size_dd=fv_fm_size_dd, xl_df=xl_df, param='yield', sp=sp)
            self._populate_fv_fm_sample_rows(ind, fv_fm_size_dd, xl_df, 'cylinder_vol', sp=sp)
        return fv_fm_size_dd

    def _populate_recruit_survival(self, sp, survival_data):
        xl_df = pd.read_excel(io=self.physiological_data_path, sheet_name=f'{sp}_recruit_survival_bh')
        self._pop_survival_data(survival_data=survival_data, xl_df=xl_df, species=sp,
                                adult_recruit='recruit', time_unit='month', exp_type='main')

    def _populate_adult_survival(self, sp, survival_data):
        xl_df = pd.read_excel(io=self.physiological_data_path, sheet_name=f'{sp}_adult_survival_bh')
        self._pop_survival_data(survival_data=survival_data, xl_df=xl_df, species=sp,
                                adult_recruit='adult', time_unit='day', exp_type='main')

    def _populate_fv_fm_sample_rows(self, ind, fv_fm_size_dd, xl_df, param, sp):
        # The rack number is not so straight forward to ascertain.
        # For the ad and ah we can use the first character alone of the 'rack' value
        # For the pd we need to extract the number from the 'rack' value.
        if sp == 'pd':
            rack_match = re.compile('\d+')
            rack = rack_match.findall(xl_df.at[ind, 'rack'])[0]
            if len(rack) > 1:
                foo = 'bar'
        elif sp in ['ad', 'ah', 'lp']:
            rack = xl_df.at[ind, 'rack'][0]

        ident = f"{xl_df.at[ind, 'temp']}_{xl_df.at[ind, 'tank']}_{rack}_" \
                f"{xl_df.at[ind, 'rack_row']}_{xl_df.at[ind, 'rack_col']}_{xl_df.at[ind, 'time']}"
        if not np.isnan(xl_df.at[ind, param]):
            # Then there is a valid yield for this sample
            if ident in fv_fm_size_dd.keys():
                # Check to see if valid param already exists
                self._check_valid_param_val_doesnt_exist(ident, param, fv_fm_size_dd, sp)
            else:
                # no need to check
                pass
            fv_fm_size_dd[ident][param] = float(xl_df.at[ind, param])

            if 'time' in fv_fm_size_dd[ident].keys():
                assert(xl_df.at[ind, 'time'] == fv_fm_size_dd[ident]['time'])
            else:
                fv_fm_size_dd[ident]['time'] = xl_df.at[ind, 'time']
        else:
            fv_fm_size_dd[ident][param] = np.nan
            if 'time' in fv_fm_size_dd[ident].keys():
                assert (xl_df.at[ind, 'time'] == fv_fm_size_dd[ident]['time'])
            else:
                fv_fm_size_dd[ident]['time'] = xl_df.at[ind, 'time']

    def _check_valid_param_val_doesnt_exist(self, ident, param, fv_fm_size_dd, sp):
        if param not in fv_fm_size_dd[ident].keys():
            return
        elif np.isnan(fv_fm_size_dd[ident][param]):
            return
        else:
            print(f'A valid {param} value already exists for {ident} for {sp}')


    def _pop_survival_data(self, survival_data, xl_df, species, adult_recruit, time_unit, exp_type):
        for row in range(len(xl_df.index)):
            for col in range(2, len(list(xl_df)), 1):
                survival_data.append([
                # species
                species,
                # adult_recruit
                adult_recruit,
                # timve_value
                list(xl_df)[col],
                # time_unit
                time_unit,
                # tank
                xl_df['tank'].iat[row],
                # temp
                xl_df['temp'].iat[row],
                # exp_type
                exp_type,
                # survival
                xl_df.iat[row, col]
                            ])

schupp_figures().plot_fig_1()

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

TODO we need to split up the contents of this figure and create some new figures.
We will refactor the original class a base class and then call the creation of the figures
as seperate clases.
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
import itertools
import sys
import random
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
from datetime import datetime
DEGREE_SIGN = u'\N{DEGREE SIGN}'
from sputils.spbars import SPBars
from sputils.sphierarchical import SPHierarchical
from scipy import stats


class SchuppFigures:
    def __init__(self):
        # Paths
        self.root_dir = os.path.abspath(os.path.dirname(__file__))
        self.sp_datasheet_path = os.path.join(self.root_dir, '20201018_schupp_all_data_datasheet.xlsx')
        self.sp_data_path = os.path.join(self.root_dir, '20201018_schupp_all')
        self.physiological_data_path = os.path.join(self.root_dir, 'schupp_et_al_physiological_coral_data.xlsx')
        self.sp_seq_count_path = os.path.join(
            self.sp_data_path,
            'post_med_seqs',
            '125_20201018_DBV_20201020T020625.seqs.absolute.abund_and_meta.txt'
        )
        self.sp_profile_abund_and_meta_path = os.path.join(
            self.sp_data_path, 'its2_type_profiles',
            '125_20201018_DBV_20201020T020625.profiles.absolute.abund_and_meta.txt'
        )
        self.sp_profile_abund_path = os.path.join(
            self.sp_data_path, 'its2_type_profiles',
            '125_20201018_DBV_20201020T020625.profiles.absolute.abund_only.txt'
        )
        self.sp_profile_meta_path = os.path.join(
            self.sp_data_path, 'its2_type_profiles',
            '125_20201018_DBV_20201020T020625.profiles.meta_only.txt'
        )
        self.sp_between_smp_dist_path_c = os.path.join(
            self.sp_data_path,
            'between_sample_distances', 'C',
            '20201020T020625_braycurtis_sample_distances_C_sqrt.dist'
        )
        self.sp_between_smp_dist_path_d = os.path.join(
            self.sp_data_path,
            'between_sample_distances', 'D',
            '20201020T020625_braycurtis_sample_distances_D_sqrt.dist'
        )

        # Datatypes and species
        self.species_full = [
            'Acropora digitifera', 'Acropora hyacinthus', 'Pocillopora damicornis', 'Leptastrea purpurea'
        ]
        self.species_short = ['ad', 'ah', 'pd', 'lp']
        self.species_short_to_full_dict = {
            'ad':'Acropora digitifera', 'ah':'Acropora hyacinthus',
            'pd':'Pocillopora damicornis', 'lp':'Leptastrea purpurea'
        }
        self.data_types = ['adult_survival', 'recruit_survival', 'recruit_size', 'recruit_fv_fm']

        # Dataframe
        # We want to create one single dataframe that holds all of the physiology data.
        # We will hold the zooxs data in seperate dataframes that we will create from the SP outputs
        # For the physiology data we will need a dataframe for survival and one for growth and fv_fm
        # This is because survival was assesed as a batch (i.e. 3 dead out of 5), while the growth and fv_fm
        # were done on a persample basis.

        fv_fm_size_data, survival_data = self._populate_data_holders()
        self.survival_df, self.fv_fm_size_df = self._make_survival_fv_fm_size_dfs(fv_fm_size_data, survival_data)
        self.sp_sample_uid_to_sample_name_dict = None
        self.sp_sample_name_to_sample_uid_dict = None
        self.sp_profile_uid_to_profile_name_dict = None
        self.sp_profile_name_to_profile_uid_dict = None
        self.sp_datasheet_df, self.sp_seq_abs_abund_df, self.sp_profile_meta_df, self.sp_profile_abs_abund_df = self._get_sp_dfs()
        # Make the relative abundance dfs
        self.sp_seq_rel_abund_df = self.sp_seq_abs_abund_df.div(
            self.sp_seq_abs_abund_df.sum(axis=1), axis=0
        )
        # Create a dictionary that has the sample_uid to relative proportion of the sequence D2d which
        # we will use for the D clustering
        self.d2d_to_sample_uid = dict(zip(self.sp_seq_rel_abund_df.index.values, self.sp_seq_rel_abund_df.D2d))
        self.sp_profile_rel_abund_df = self.sp_profile_abs_abund_df.div(
            self.sp_profile_abs_abund_df.sum(axis=1), axis=0
        )

        # Get a useful list of the SymPortal sample names and sample UIDs that are used
        # in the study.
        # We use this because there are samples in the SymPortal submission that are not used in this study
        # So it is good to know which samples we need to be working with
        self.sp_sample_names_of_study, self.sp_sample_uids_of_study = self._get_name_and_uid_of_samples_to_plot()

        # Reorder the abundance dataframes in order of most abundant sequences/profiles
        sorted_seq_index = self.sp_seq_rel_abund_df.sum(axis=0).sort_values(ascending=False).index
        sorted_profile_index = self.sp_profile_rel_abund_df.sum(axis=0).sort_values(ascending=False).index
        self.sp_seq_rel_abund_df = self.sp_seq_rel_abund_df.reindex(sorted_seq_index, axis=1)
        self.sp_seq_abs_abund_df = self.sp_seq_abs_abund_df.reindex(sorted_seq_index, axis=1)
        self.sp_profile_rel_abund_df = self.sp_profile_rel_abund_df.reindex(sorted_profile_index, axis=1)
        self.sp_profile_abs_abund_df = self.sp_profile_abs_abund_df.reindex(sorted_profile_index, axis=1)

    def _get_name_and_uid_of_samples_to_plot(self):
        sample_names_list = []
        # recruit ad
        sample_names_list += list(self.sp_datasheet_df[
                                      (self.sp_datasheet_df['host_species'] == 'digitifera') &
                                      (self.sp_datasheet_df['age'].str.contains("month")) &
                                      (~self.sp_datasheet_df['age'].str.contains("through")) &
                                      (~self.sp_datasheet_df['age'].str.contains("delayed"))
                                      ].index.values)
        # recruit ah
        sample_names_list += list(self.sp_datasheet_df[
                                      (self.sp_datasheet_df['host_species'] == 'surculosa') &
                                      (self.sp_datasheet_df['age'].str.contains("month")) &
                                      (~self.sp_datasheet_df['age'].str.contains("adult"))
                                      ].index.values)
        # recruit pd
        sample_names_list += list(self.sp_datasheet_df[
                                      (self.sp_datasheet_df['host_species'] == 'damicornis') &
                                      (self.sp_datasheet_df['age'].str.contains("month")) &
                                      (~self.sp_datasheet_df['age'].str.contains("adult"))
                                      ].index.values)
        # recruit lp
        sample_names_list += list(self.sp_datasheet_df[
                                      (self.sp_datasheet_df['host_species'] == 'purpurea') &
                                      (self.sp_datasheet_df['age'].str.contains("month")) &
                                      (~self.sp_datasheet_df['age'].str.contains("adult"))
                                      ].index.values)
        # adult ad
        sample_names_list += list(self.sp_datasheet_df[
                                      (self.sp_datasheet_df['age'] == 'adult') &
                                      (self.sp_datasheet_df['host_species'] == 'digitifera')
                                      ].index.values)
        # adult ah
        sample_names_list += list(self.sp_datasheet_df[
                                      (self.sp_datasheet_df['age'] == 'adult') &
                                      (self.sp_datasheet_df['host_species'] == 'surculosa')
                                      ].index.values)
        # adult pd
        sample_names_list += list(self.sp_datasheet_df[
                                      (self.sp_datasheet_df['age'] == 'adult') &
                                      (self.sp_datasheet_df['host_species'] == 'damicornis')
                                      ].index.values)
        # adult lp
        sample_names_list += [_ for _ in self.sp_datasheet_df[
            (self.sp_datasheet_df['age'] == 'adult') &
            (self.sp_datasheet_df['host_species'] == 'purpurea')
            ].index.values if
                              "P4" in _]

        return sample_names_list, [self.sp_sample_name_to_sample_uid_dict[_] for _ in sample_names_list]

    def _get_sp_dfs(self):
        """
        Return three pd DataFrames.
        One from the SP datasheet
        One that contains the post-MED seq counts
        One that contains the profile counts
        """
        sp_ds_df = self._make_sp_datasheet_df()

        seq_count_df = self._make_sp_seq_abund_df()

        profile_abund_df, profile_meta_df = self._make_profile_dfs()

        return sp_ds_df, seq_count_df, profile_meta_df, profile_abund_df

    def _make_profile_dfs(self):
        profile_meta_df = self._make_sp_profile_meta_df()
        profile_abund_df = self._make_sp_profile_abund_df()
        return profile_abund_df, profile_meta_df

    def _make_sp_profile_abund_df(self):
        profile_abund_df = pd.read_csv(self.sp_profile_abund_path, sep='\t')
        profile_abund_df.set_index("sample_uid", drop=True, inplace=True)
        return profile_abund_df

    def _make_sp_profile_meta_df(self):
        profile_meta_df = pd.read_csv(self.sp_profile_meta_path, sep='\t')
        self.sp_profile_uid_to_profile_name_dict = dict(zip(
            profile_meta_df["ITS2 type profile UID"],
            profile_meta_df['ITS2 type profile']
        ))
        self.sp_profile_name_to_profile_uid_dict = dict(zip(
            profile_meta_df['ITS2 type profile'],
            profile_meta_df["ITS2 type profile UID"]
        ))
        profile_meta_df.rename(
            columns={"ITS2 type profile": "profile_name", "ITS2 type profile UID": "profile_uid"},
            inplace=True
        )
        profile_meta_df.set_index("profile_uid", inplace=True, drop=True)
        return profile_meta_df

    def _make_sp_seq_abund_df(self):
        seq_count_df = pd.read_csv(self.sp_seq_count_path, sep='\t')
        seq_count_df = seq_count_df.iloc[:-1, :]
        seq_count_df.sample_uid = seq_count_df.sample_uid.astype(int)
        seq_count_df.set_index(keys='sample_uid', drop=True, inplace=True)
        self.sp_sample_uid_to_sample_name_dict = dict(zip(seq_count_df.index.values, seq_count_df.sample_name))
        self.sp_sample_name_to_sample_uid_dict = dict(zip(seq_count_df.sample_name, seq_count_df.index.values))
        index_of_first_seq = list(seq_count_df).index("A3")
        seq_count_df = seq_count_df.iloc[:-1, index_of_first_seq:]
        self.c_sample_uid_to_post_med_absolute_dict, self.c_sample_uid_to_post_med_unique_dict = self._make_clade_specific_post_med_absolute_and_unique_dicts(clade='C', seq_count_df=seq_count_df)
        self.d_sample_uid_to_post_med_absolute_dict, self.d_sample_uid_to_post_med_unique_dict = self._make_clade_specific_post_med_absolute_and_unique_dicts(
            clade='D', seq_count_df=seq_count_df)
        return seq_count_df

    def _make_clade_specific_post_med_absolute_and_unique_dicts(self, clade, seq_count_df):
        abs_dict = {}
        unique_dict = {}
        for sample_ind in seq_count_df.index:
            ser = seq_count_df.loc[sample_ind][list([_ for _ in seq_count_df if (_.startswith(clade) or _.endswith(clade))])]
            ser = ser[ser != 0]
            abs_dict[sample_ind] = ser.sum()
            unique_dict[sample_ind] = len(ser)
        return abs_dict, unique_dict


    def _make_sp_datasheet_df(self):
        sp_ds_df = pd.read_excel(io=self.sp_datasheet_path, skiprows=[0], engine='openpyxl')
        collection_cols = [lab for lab in list(sp_ds_df) if ("collection" in lab) or ("Unnamed" in lab)]
        sp_ds_df.drop(columns=collection_cols, inplace=True)
        sp_ds_df.set_index("sample_name", drop=True, inplace=True)
        return sp_ds_df

    def _make_survival_fv_fm_size_dfs(self, fv_fm_size_data, survival_data):
        survival_df = pd.DataFrame(
            data=survival_data,
            columns=['species', 'adult_recruit', 'time_value', 'time_unit', 'tank', 'temperature', 'exp_type',
                     'survival'])
        fv_fm_size_df = pd.DataFrame(
            data=fv_fm_size_data,
            columns=[
                'species', 'adult_recruit', 'time_value', 'time_unit', 'temperature', 'tank', 'rack',
                'rack_row', 'rack_col', 'cyl_vol', 'fv_fm', 'exp_type',
            ])

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

                elif data_type == 'recruit_size':
                    fv_fm_size_dd = self._populate_recruit_fv_fm_size_dd(sp)
                    self._populate_recruit_fv_fm_size_data(fv_fm_size_data, fv_fm_size_dd, sp)
        return fv_fm_size_data, survival_data

    def _populate_recruit_fv_fm_size_data(self, fv_fm_size_data, fv_fm_size_dd, sp):
        # At this point we have the dict populated and we can now populate the fv_fm_size_data
        for k in fv_fm_size_dd.keys():
            temp, tank, rack, rack_row, rack_col, time = k.split('_')
            fv_fm_size_data.append([
                sp, 'recruit', int(time), 'months',
                temp, tank, rack, rack_row, rack_col,
                fv_fm_size_dd[k]['cylinder_vol'], fv_fm_size_dd[k]['yield'], 'main'
            ])

    def _populate_recruit_fv_fm_size_dd(self, sp):
        # A dictionary where identifier is key and value is a dict with keys of 'fv_fm' and 'cyl_vol
        # Identifier will be made in the form "{temp}_{tank}_{rack}_{rack_row}_{rack_row}"
        # We will place numpy.np for val if no val is available
        # We will check to see that there are not duplicate entries for a given sample and data_type
        fv_fm_size_dd = defaultdict(dict)
        xl_df = pd.read_excel(io=self.physiological_data_path, sheet_name=f'{sp}_recruit_size_fv_fm_bh', engine='openpyxl')
        for ind in xl_df.index:
            self._populate_fv_fm_sample_rows(ind=ind, fv_fm_size_dd=fv_fm_size_dd, xl_df=xl_df, param='yield', sp=sp)
            self._populate_fv_fm_sample_rows(ind, fv_fm_size_dd, xl_df, 'cylinder_vol', sp=sp)
        return fv_fm_size_dd

    def _populate_recruit_survival(self, sp, survival_data):
        xl_df = pd.read_excel(io=self.physiological_data_path, sheet_name=f'{sp}_recruit_survival_bh', engine='openpyxl')
        self._pop_survival_data(survival_data=survival_data, xl_df=xl_df, species=sp,
                                adult_recruit='recruit', time_unit='months', exp_type='main')

    def _populate_adult_survival(self, sp, survival_data):
        xl_df = pd.read_excel(io=self.physiological_data_path, sheet_name=f'{sp}_adult_survival_bh', engine='openpyxl')
        self._pop_survival_data(survival_data=survival_data, xl_df=xl_df, species=sp,
                                adult_recruit='adult', time_unit='days', exp_type='main')

    def _populate_fv_fm_sample_rows(self, ind, fv_fm_size_dd, xl_df, param, sp):
        # The rack number is not so straight forward to ascertain.
        # For the ad and ah we can use the first character alone of the 'rack' value
        # For the pd we need to extract the number from the 'rack' value.
        ident = self._get_ident(ind, sp, xl_df)
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
            # Check to see if param value already exists
            if param in fv_fm_size_dd[ident].keys():
                # Already an entry present
                pass
            else:
                # No entry present
                fv_fm_size_dd[ident][param] = np.nan
            if 'time' in fv_fm_size_dd[ident].keys():
                assert (xl_df.at[ind, 'time'] == fv_fm_size_dd[ident]['time'])
            else:
                fv_fm_size_dd[ident]['time'] = xl_df.at[ind, 'time']

    def _get_ident(self, ind, sp, xl_df):
        if sp == 'pd':
            rack_match = re.compile('\d+')
            rack = rack_match.findall(xl_df.at[ind, 'rack'])[0]
            if len(rack) > 1:
                foo = 'bar'
        elif sp in ['ad', 'ah', 'lp']:
            rack = xl_df.at[ind, 'rack'][0]
        ident = f"{xl_df.at[ind, 'temp']}_{xl_df.at[ind, 'tank']}_{rack}_" \
                f"{xl_df.at[ind, 'rack_row']}_{xl_df.at[ind, 'rack_col']}_{xl_df.at[ind, 'time']}"
        return ident

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

class PhysiologicalBasePlot(SchuppFigures):
    """
    Plot a figure that has four columns one for each specie of coral
    and 4 rows one for each data type:
    adult survival
    recruit survival
    size
    fv/fm
    """
    def __init__(self):
        super().__init__()
        # Figure
        # 10 wide, 6 deep
        self.fig = plt.figure(figsize=(10, 10))
        # 6 down 4 across
        self.gs = gridspec.GridSpec(4, 4)
        self.axes = [[] for sp in self.species_short]

        for i, sp in enumerate(self.species_short):
            for j, data_type in enumerate(self.data_types):
                self.axes[i].append(plt.subplot(self.gs[j, i]))

        self.temp_plot_marker_dict = {29: "v", 30: "o", 31: "^", '29': "v", '30': "o", '31': "^"}
        self.temp_plot_colour_dict = {29: "#a6a6a6", 30: "#696969", 31: "#0d0d0d", '29': "#a6a6a6", '30': "#696969",
                                      '31': "#0d0d0d"}
        self.marker_size = 4

    def plot(self):
        """
        The main plotting method for plotting up the figure that has the four species columns with a plot
        for each of the data types.
        """
        for i, sp in enumerate(self.species_short):
            for j, data_type in enumerate(self.data_types):
                if data_type == 'adult_survival':
                    self._plot_adult_survival(i, j, sp, data_type)
                elif data_type == 'recruit_survival':
                    self._plot_recruit_survival(data_type, i, j, sp)
                elif data_type == 'recruit_size':
                    self._plot_recruit_size(data_type, i, j, sp)
                elif data_type == 'recruit_fv_fm':
                    self._plot_recruit_fv_fm(data_type, i, j, sp)
        print('saving svg')
        plt.savefig(os.path.join(self.root_dir, 'figures',
                                 f"physiology_and_survival_base_figure_{str(datetime.now()).split('.')[0].replace('-', '').replace(' ', 'T').replace(':', '')}.svg"),
                    dpi=1200)
        print('saving png')
        plt.savefig(os.path.join(self.root_dir, 'figures',
                                 f"physiology_and_survival_base_figure_{str(datetime.now()).split('.')[0].replace('-', '').replace(' ', 'T').replace(':', '')}.png"),
                    dpi=1200)

    def _plot_recruit_fv_fm(self, data_type, i, j, sp):
        ax = self.axes[i][j]
        param_to_plot = "fv_fm"
        working_df = self.fv_fm_size_df[
            (self.fv_fm_size_df['species'] == sp) & (self.fv_fm_size_df[param_to_plot].notnull())]
        self._plot_if_working_df(ax, data_type, param_to_plot, sp, working_df)

    def _plot_if_working_df(self, ax, data_type, param_to_plot, sp, working_df, time_first_after_val=None):
        if working_df.index.to_list():
            self._plot_a_set_of_line_data(ax, data_type, sp, working_df, param_to_plot, time_first_after_val)
        else:
            # Turn off the axis but make sure that a title is still provided if this is adult survival
            ax.axis('off')
            if data_type == 'adult_survival':
                ax.set_title(self.species_short_to_full_dict[sp], fontsize='small', style='italic')


    def _plot_recruit_size(self, data_type, i, j, sp):
        ax = self.axes[i][j]
        param_to_plot = "cyl_vol"
        working_df = self.fv_fm_size_df[
            (self.fv_fm_size_df['species'] == sp) & (self.fv_fm_size_df[param_to_plot].notnull())].copy()
        # Convert to cm3
        working_df['cyl_vol'] = working_df['cyl_vol'] / 1000
        self._plot_if_working_df(ax, data_type, param_to_plot, sp, working_df)

    def _plot_recruit_survival(self, data_type, i, j, sp):
        # This is slightly more complicated because we need to work with the shipping of the corals from Guam to Germany
        # To highlight this I think it will be a good idea to have a vertical
        # line or to shade the backgrounds a different colour
        ax = self.axes[i][j]
        # We will need to do an ax.errorbar for each of the temperatures
        param_to_plot = "survival_percent"
        working_df = self.survival_df[
            (self.survival_df['adult_recruit'] == 'recruit') &
            (self.survival_df['species'] == sp)
        ].copy()
        # The survival percent calculation will depend on the species because there were slightly
        # different time points for each of the species.
        # For the time points measured before moving, the survival will be the value divided by
        # the starting value * 100
        # For ad, ah and lp any time value < 4 is a before shipment point and any => 4 is after.
        # For pd <3 is before and >= 3 is after

        # first get the starting survivals and plug these into a dict
        # There will be one per tank
        # Pull out the time 0s into a series
        # The for each row get a tank
        # Then make these k, v pairs in dict
        time_zero_dict = self._make_time_zero_df_recruit_survival(working_df)

        # First create a dict that holds the survial percentages for the last measurement before
        # the move
        time_first_after_val, time_last_before_dict = self._make_time_last_before_dict_recruit_survival(
            sp, time_zero_dict, working_df)

        # Secondly, we want a dict that given a tank we will be able to get the
        # time_first_after_survival. I.e. the first survival point after the move
        time_first_after_survival_dict = self._make_time_first_after_survival_dict(time_first_after_val,
                                                                                   working_df)

        # Now that we have the time_last_before_dict we will be able to work out the survival_percentage
        # for the samples after the move

        # Finally calculate the survival percentages for all samples
        self._calculate_survival_percent_recruit_survival(time_first_after_survival_dict,
                                                          time_first_after_val, time_last_before_dict,
                                                          time_zero_dict, working_df)

        working_df = working_df[working_df[param_to_plot].notnull()]
        # Now we can finally plot
        # Similar to the adults we want to do an ax.errorbar for each of the temperatures
        # NB for the recruit survival points, we don't want there to be a line connecting the before
        # and after travel points.
        self._plot_if_working_df(ax, data_type, param_to_plot, sp, working_df, time_first_after_val)

    def _calculate_survival_percent_recruit_survival(self, time_first_after_survival_dict, time_first_after_val,
                                                     time_last_before_dict, time_zero_dict, working_df):
        survival_percent = []
        for ind in working_df.index:
            tank = working_df.at[ind, 'tank']
            if working_df.at[ind, 'time_value'] < time_first_after_val:
                # Then we want to calculate against the start value
                survival_percent.append((working_df.at[ind, 'survival'] / time_zero_dict[tank]) * 100)
            else:
                # The percentage survial for the given sample at the last time point before move
                time_last_before_survival = time_last_before_dict[tank]
                # The absolute survial value for the given sample at the first time point after move
                time_first_after_survival = time_first_after_survival_dict[tank]
                # The decrease in survial for current time point from the last time point before move
                survival_decrease_percent_after = ((time_first_after_survival - working_df.at[
                    ind, 'survival']) / time_first_after_survival) * 100
                # Finally, the survival percentage is a little complicated. Hard to explain. Just think about it.
                survival_percent.append(time_last_before_survival - (time_last_before_survival*(survival_decrease_percent_after/100)))
        working_df['survival_percent'] = survival_percent

    def _make_time_first_after_survival_dict(self, time_first_after_val, working_df):
        time_first_after_survival_dict = {}
        time_first_after_df = working_df[working_df['time_value'] == time_first_after_val]
        for ind in time_first_after_df.index:
            tank = time_first_after_df.at[ind, 'tank']
            time_first_after_survival_dict[tank] = time_first_after_df.at[ind, 'survival']
        return time_first_after_survival_dict

    def _make_time_last_before_dict_recruit_survival(self, sp, time_zero_dict, working_df):
        time_last_before_dict = {}
        if sp in ['ad', 'ah']:
            time_last_before_val = 3
            time_first_after_val = 4
        elif sp == 'lp':
            time_last_before_val = 3.9
            time_first_after_val = 4
        elif sp == 'pd':
            time_last_before_val = 2.9
            time_first_after_val = 3
        else:
            raise RuntimeError(f'Unexpected species short {sp}')
        time_last_before_df = working_df[working_df['time_value'] == time_last_before_val]
        for ind in time_last_before_df.index:
            tank = time_last_before_df.at[ind, 'tank']
            time_last_before_dict[tank] = (time_last_before_df.at[ind, 'survival'] / time_zero_dict[tank]) * 100
        return time_first_after_val, time_last_before_dict

    def _make_time_zero_df_recruit_survival(self, working_df):
        time_zero_df = working_df[working_df['time_value'] == 0]
        time_zero_dict = {}
        for ind in time_zero_df.index:
            tank = self.survival_df.at[ind, 'tank']
            time_zero_dict[tank] = time_zero_df.at[ind, 'survival']
        return time_zero_dict

    def _plot_a_set_of_line_data(
            self, ax, data_type, sp, working_df, param_to_plot, time_first_after_val=None
    ):
        for temp in working_df['temperature'].unique():
            if data_type == 'recruit_survival':
                # Plot the points in two seperate calls, one for before shipping and one for after so that
                # there is a break in the line between the points
                ser = working_df[(working_df['temperature'] == temp) & (working_df['time_value'] < time_first_after_val)]
                # Calc average survival for each time point and the standard error of the mean
                self._calc_mean_sem_plot_line(ax, param_to_plot, ser, temp)
                ser = working_df[(working_df['temperature'] == temp) & (working_df['time_value'] >= time_first_after_val)]
                # Calc average survival for each time point and the standard error of the mean
                self._calc_mean_sem_plot_line_no_label(ax, param_to_plot, ser, temp)
            else:
                ser = working_df[working_df['temperature'] == temp]
                # Calc average survival for each time point and the standard error of the mean
                self._calc_mean_sem_plot_line(ax, param_to_plot, ser, temp)
        if "survival" in data_type:
            ax.set_ylim(-10,110)
        time_unit = working_df['time_unit'].unique()[0]
        if time_unit == 'days':
            ax.set_xticks([0,5,10,15,20,25])
        elif time_unit == 'months':
            ax.set_xticks([0,3,6,9,12])
        if sp == 'ad' and data_type == 'adult_survival':
            ax.legend(loc='lower left', fontsize='xx-small')
            ax.set_ylabel('Adult survial %', fontsize='small')
        ax.set_xlabel(time_unit, fontsize='xx-small')
        if sp == 'ad' and data_type == 'recruit_survival':
            # only set one per row
            ax.set_ylabel('Recruit survial %', fontsize='small')
        elif sp == 'ad' and "size" in data_type:
            ax.set_ylabel('Size\ncyl. vol. ml', fontsize='small')
        elif sp == 'ad' and data_type == "recruit_fv_fm":
            ax.set_ylabel('Fv/Fm', fontsize='small')
        if data_type == 'adult_survival':
            ax.set_title(self.species_short_to_full_dict[sp], fontsize='small', style='italic')
        if data_type == "recruit_fv_fm":
            ax.set_ylim(0.54,0.71)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.tight_layout()

    def _calc_mean_sem_plot_line(self, ax, param_to_plot, ser, temp):
        means = [ser[ser['time_value'] == time_val][param_to_plot].mean() for time_val in
                 ser['time_value'].unique()]
        sem = [ser[ser['time_value'] == time_val][param_to_plot].sem() for time_val in
               ser['time_value'].unique()]
        ax.errorbar(
            x=ser['time_value'].unique(), y=means, yerr=sem, marker=self.temp_plot_marker_dict[temp],
            linestyle='--', linewidth=1, ecolor=self.temp_plot_colour_dict[temp],
            elinewidth=1, color=self.temp_plot_colour_dict[temp], markersize=self.marker_size,
            label=f'{temp}{DEGREE_SIGN}C'
        )

    def _calc_mean_sem_plot_line_no_label(self, ax, param_to_plot, ser, temp):
        means = [ser[ser['time_value'] == time_val][param_to_plot].mean() for time_val in
                 ser['time_value'].unique()]
        sem = [ser[ser['time_value'] == time_val][param_to_plot].sem() for time_val in
               ser['time_value'].unique()]
        ax.errorbar(
            x=ser['time_value'].unique(), y=means, yerr=sem, marker=self.temp_plot_marker_dict[temp],
            linestyle='--', linewidth=1, ecolor=self.temp_plot_colour_dict[temp],
            elinewidth=1, color=self.temp_plot_colour_dict[temp], markersize=self.marker_size,
            label=None
        )

    def _plot_adult_survival(self, i, j, sp, data_type):

        # Plot the aduult_survival as just that, survival. So start with 100% and then decrease
        ax = self.axes[i][j]
        # We will need to do an ax.errorbar for each of the temperatures
        param_to_plot = "survival_percent"
        working_df = self.survival_df[
            (self.survival_df['adult_recruit'] == 'adult') & (self.survival_df['species'] == sp)].copy()
        working_df['survival_percent'] = (working_df['survival'] / 5) * 100
        working_df = working_df[working_df[param_to_plot].notnull()]
        self._plot_if_working_df(ax, data_type, param_to_plot, sp, working_df)

class HierarchicalPlot(SchuppFigures):
    """
    Plot that will show the two clusters of the between sample distances
    One plot for Cladocopium and one for Durusdinium
    Each plot will be made of 3 plots
    Top will be the hierarchical cluster plot
    Middle will be the sequences and profiles in order of the hierarchical plot
    Bottom will indicate the artiical clustering by us
    bottom again will be host species

    We will also produce a seperate set of supporting figure that are
    histograms of the the post_med_absolute and post_med_unique (clade separated)
    sequences for the samples. We will use these to justify the filtering
    of samples that we will implement for the main figure.
    For D will not use those samples with <10 unique D seqs and <5000 absolute sequences
    """
    def __init__(self):
        super().__init__()
        # Get the D samples to plot
        self.d_sph_no_plot = SPHierarchical(dist_output_path=self.sp_between_smp_dist_path_d, no_plotting=True)
        self.d_sample_uids_in_dist = list(self.d_sph_no_plot.dist_df)
        # we want to be working with the following intersect
        self.d_sample_uids_to_plot_non_filtered = [_ for _ in self.d_sample_uids_in_dist if _ in self.sp_sample_uids_of_study]

        # Get the D samples to plot
        self.c_sph_no_plot = SPHierarchical(dist_output_path=self.sp_between_smp_dist_path_c, no_plotting=True)
        self.c_sample_uids_in_dist = list(self.c_sph_no_plot.dist_df)
        # we want to be working with the following intersect
        self.c_sample_uids_to_plot_non_filtered = [_ for _ in self.c_sample_uids_in_dist if
                                                   _ in self.sp_sample_uids_of_study]

    def plot_supporting_histograms(self):
        fig = plt.figure(figsize=(10, 10))
        # 6 down 4 across
        gs = gridspec.GridSpec(2, 2)
        axes = []
        plot_tyes = ['post_med_absolute', 'post_med_unique']
        for i, genus in enumerate(['C', 'D']):
            for j, plot_type in enumerate(plot_tyes):
                axes.append(plt.subplot(gs[i:i+1, j:j+1]))

        # D absolute
        d_post_med_absolute_values = [self.d_sample_uid_to_post_med_absolute_dict[_] for _ in self.d_sample_uids_to_plot_non_filtered]
        d_abs_kde = stats.gaussian_kde(d_post_med_absolute_values)
        d_abs_bins = range(0, 50000, int(50000 / 20))
        d_abs_kde_x = np.linspace(0, 50000, 100)
        axes[0].hist(d_post_med_absolute_values, bins=d_abs_bins, density=True)
        axes[0].plot(d_abs_kde_x, d_abs_kde(d_abs_kde_x))
        axes[0].set_title('Durusdinium: post-MED absolute sequences')
        axes[0].set_xlabel('post-MED absolute sequences')

        # D unique
        d_post_med_unique_values = [self.d_sample_uid_to_post_med_unique_dict[_] for _ in self.d_sample_uids_to_plot_non_filtered]
        d_unique_kde = stats.gaussian_kde(d_post_med_unique_values)
        d_unique_bins = range(0, 40, 2)
        d_unique_kde_x = np.linspace(0, 40, 100)
        axes[1].hist(d_post_med_unique_values, bins=d_unique_bins, density=True)
        axes[1].plot(d_unique_kde_x, d_unique_kde(d_unique_kde_x))
        axes[1].set_title('Durusdinium: post-MED unique sequences')
        axes[1].set_xlabel('post-MED unique sequences')

        # C absolute
        c_post_med_absolute_values = [self.c_sample_uid_to_post_med_absolute_dict[_] for _ in
                                        self.c_sample_uids_to_plot_non_filtered]
        c_abs_kde = stats.gaussian_kde(c_post_med_absolute_values)
        c_abs_bins = range(0, 50000, int(50000 / 20))
        c_abs_kde_x = np.linspace(0, 50000, 100)
        axes[2].hist(c_post_med_absolute_values, bins=c_abs_bins, density=True)
        axes[2].plot(c_abs_kde_x, c_abs_kde(c_abs_kde_x))
        axes[2].set_title('Cladocopium: post-MED absolute sequences')
        axes[2].set_xlabel('post-MED absolute sequences')

        # C unique
        c_post_med_unique_values = [self.c_sample_uid_to_post_med_unique_dict[_] for _ in
                                  self.c_sample_uids_to_plot_non_filtered]
        c_unique_kde = stats.gaussian_kde(c_post_med_unique_values)
        c_unique_bins = range(0, 40, 2)
        c_unique_kde_x = np.linspace(0, 40, 100)
        axes[3].hist(c_post_med_unique_values, bins=c_unique_bins, density=True)
        axes[3].plot(c_unique_kde_x, c_unique_kde(c_unique_kde_x))
        axes[3].set_title('Cladocopium: post-MED unique sequences')
        axes[3].set_xlabel('post-MED unique sequences')

        plt.savefig(os.path.join(self.root_dir, 'figures',
                                 f"sup_histograms_{str(datetime.now()).split('.')[0].replace('-', '').replace(' ', 'T').replace(':', '')}.svg"),
                    dpi=1200)
        plt.savefig(os.path.join(self.root_dir, 'figures',
                                 f"sup_histograms_{str(datetime.now()).split('.')[0].replace('-', '').replace(' ', 'T').replace(':', '')}.png"),
                    dpi=1200)

    def plot_main_hierarchical_clutering_figure(self):
        # Get the list of samples that need plotting
        # self.sp_between_smp_dist_path_c
        # self.sp_between_smp_dist_path_d
        self.fig = plt.figure(figsize=(10, 10))
        # 6 down 4 across
        # TODO we will need to adjust this as we refine the figure
        self.gs = gridspec.GridSpec(8, 1)
        self.axes = []
        plot_tyes = ['hierarchical', 'seq_prof', 'cluster', 'species']
        for i, genus in enumerate(['Durusdinium', 'Cladocopium']):
            for j, plot_type in enumerate(plot_tyes):
                self.axes.append(plt.subplot(self.gs[((i * len(plot_tyes)) + j):((i * len(plot_tyes)) + j) + 1, :]))
        # TODO we need to plot up the unique and aboslute post-meds on a clade split basis
        # particularly for C I expect this to get rid of a large proportion of the sampls.

        # For D it looks like we can screen out most of the bad distances using a conservative
        # filtering based on post_med_absolute and post_med_unique
        # We will justify the quantificaiton of the cutoff using a historgram and kde that can
        # be put into the supplementary materials but this should certainly be output as a figure


        # TODO I think that there are extra samples in the symportal data that we don't need to plot
        # there are 421 samples in the d dist matrix
        # We already worked out which samples we will be plotting the zooxs data for in the old
        # code so we should use a modificaiton of that code to get the list of samples
        # that we should be working with.
        # The sample names and sample uids that are used in the studies are held in these two objecsts
        # self.sp_sample_names_of_study, self.sp_sample_uids_of_study

        # TODO there are a number of samples that appear to be low quality and I want to potentially filter
        # these out using either the post_med_absolute (maybe <200) and or the post_med_unique values
        # As such I quickly want to visualise these values in the same order as the hierarchical clustering
        # TO do this I will create dictionaries of these values.
        # self.sample_uid_to_post_med_absolute and self.sample_uid_to_post_med_unique
        # TODO the d clustering it may be easiest to cluster by the presence of the D2d sequence
        axes = [*self.axes[:4]]
        dist_output_path = self.sp_between_smp_dist_path_d
        clade_list = ['D']
        screening_dict = self.d2d_to_sample_uid
        dendrogram_sample_uid_order = self._plot_for_clade(axes, clade_list, dist_output_path, screening_dict)

        # axes = [*self.axes[4:]]
        # dist_output_path = self.sp_between_smp_dist_path_c
        # clade_list = ['C']
        # screening_dict = self.d2d_to_sample_uid
        # self._plot_for_clade(axes, clade_list, dist_output_path, screening_dict)


        post_med_absolute_values = [self.d_sample_uid_to_post_med_absolute_dict[_] for _ in dendrogram_sample_uid_order]
        post_med_unique_values = [self.d_sample_uid_to_post_med_unique_dict[_] for _ in dendrogram_sample_uid_order]
        abs_kde = stats.gaussian_kde(post_med_absolute_values)
        unique_kde = stats.gaussian_kde(post_med_unique_values)
        abs_bins = range(0,50000,int(50000/20))
        unique_bins = range(0, 40, 2)
        abs_kde_x = np.linspace(0, 50000, 100)
        unique_kde_x = np.linspace(0,40,100)
        self.axes[4].hist(post_med_absolute_values, bins=abs_bins, density=True)
        self.axes[4].plot(abs_kde_x, abs_kde(abs_kde_x))
        self.axes[5].hist(post_med_unique_values, bins=unique_bins, density=True)
        self.axes[5].plot(unique_kde_x, unique_kde(unique_kde_x))


        self.axes[4].imshow(np.array(post_med_absolute_values)[np.newaxis, :], cmap="plasma", aspect="auto")
        self.axes[5].imshow(np.array(post_med_unique_values)[np.newaxis, :], cmap="plasma", aspect="auto")
        ex_absolute = [1000 if _ >= 5000 else 0 for _ in post_med_absolute_values]
        ex_unique = [1000 if _ >= 10 else 0 for _ in post_med_unique_values]
        self.axes[6].imshow(np.array(ex_absolute)[np.newaxis, :], cmap="plasma", aspect="auto")
        self.axes[7].imshow(np.array(ex_unique)[np.newaxis, :], cmap="plasma", aspect="auto")
        plt.show()
        foo = 'bar'

    def _plot_for_clade(self, axes, clade_list, dist_output_path, screening_dict):
        if 'C' in clade_list:
            fo = 'bar'
        sph_no_plot = SPHierarchical(dist_output_path=dist_output_path, no_plotting=True)
        sample_uids_in_dist = list(sph_no_plot.dist_df)
        # we want to be working with the following intersect
        sample_uids_to_plot = [_ for _ in sample_uids_in_dist if _ in self.sp_sample_uids_of_study]
        sample_names_to_plot = [self.sp_sample_uid_to_sample_name_dict[_] for _ in sample_uids_to_plot]
        sph_plot = SPHierarchical(dist_output_path=dist_output_path, ax=axes[0],
                                  sample_uids_included=sample_uids_to_plot)
        sph_plot.plot()
        # Thin out the lines
        axes[0].collections[0].set_linewidth(0.5)
        dendrogram_sample_uid_order = sph_plot.dendrogram['ivl']
        # Now plot up the sequences and profiles
        spb = SPBars(
            seq_count_table_path=self.sp_seq_count_path,
            profile_count_table_path=self.sp_profile_abund_and_meta_path,
            plot_type='seq_and_profile', orientation='h', legend=False, relative_abundance=True,
            sample_uids_included=dendrogram_sample_uid_order, bar_ax=axes[1], limit_genera=clade_list,
            seq_profile_scalar=(1.0, 0.3)
        )
        spb.plot()
        axes[1].set_xticks([])
        axes[1].set_yticks([])
        cluster_vals = [1 if screening_dict[_] >= 0.01 else 0 for _ in dendrogram_sample_uid_order]
        axes[2].imshow(np.array(cluster_vals)[np.newaxis, :], cmap="plasma", aspect="auto")
        species_vals = []
        for sample_name in [self.sp_sample_uid_to_sample_name_dict[_] for _ in dendrogram_sample_uid_order]:
            host_species = self.sp_datasheet_df.at[sample_name, 'host_species']
            if host_species == 'digitifera':
                species_vals.append(0)
            elif host_species == 'surculosa':
                species_vals.append(0.25)
            elif host_species == 'purpurea':
                species_vals.append(0.5)
            elif host_species == 'damicornis':
                species_vals.append(1)
        axes[3].imshow(np.array(species_vals)[np.newaxis, :], cmap="plasma", aspect="auto")
        return dendrogram_sample_uid_order

HierarchicalPlot().plot_supporting_histograms()

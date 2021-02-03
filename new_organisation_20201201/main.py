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
from sklearn.preprocessing import MinMaxScaler
import re
import itertools
import sys
import random
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
from matplotlib import cm
from datetime import datetime
DEGREE_SIGN = u'\N{DEGREE SIGN}'
from sputils.spbars import SPBars
from sputils.sphierarchical import SPHierarchical
from scipy import stats
# from matplotlib import rc
# # activate latex text rendering
# # https://stackoverflow.com/questions/8376335/styling-part-of-label-in-legend-in-matplotlib
# rc('text', usetex=True)


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
        self.sample_uid_to_d2d_rel_abund_dict = dict(zip(self.sp_seq_rel_abund_df.index.values, self.sp_seq_rel_abund_df.D2d))
        self.d_clustering_dict = {(uid): ('D1/D2d' if rel_abund >= 0.01 else 'D1/D4') for uid, rel_abund in
                             self.sample_uid_to_d2d_rel_abund_dict.items()}
        self.c_clustering_dict = {}
        for sample_ind in self.sp_seq_rel_abund_df.index:
            # if C50c makes up more than 5% of the C seqs then class as 'C50c'
            c_seqs = self.sp_seq_rel_abund_df.loc[sample_ind][[_ for _ in list(self.sp_seq_rel_abund_df) if (_.startswith('C') or _.endswith('C'))]]
            if c_seqs.sum() > 0:
                tot_c = c_seqs.sum()
                if c_seqs.C50c / tot_c > 0.15:
                    self.c_clustering_dict[sample_ind] = 'C50c'
                elif c_seqs.C66 / tot_c > 0.15:
                    self.c_clustering_dict[sample_ind] = 'C66'
                elif c_seqs.C1 / tot_c > 0.05:
                    self.c_clustering_dict[sample_ind] = 'C1'
                else:
                    self.c_clustering_dict[sample_ind] = 'other'

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

        # General colour maps
        self.species_c_map = {'digitifera':'#F8F8F8', 'surculosa': '#D0CFD4', 'damicornis': '#89888D', 'purpurea': '#4A4A4C',
                              }
        self.d_cluster_c_map = {'D1/D2d': '#104E8B', 'D1/D4': '#60AFFE', 'other': '#F8F8F8'}
        self.c_cluster_c_map = {'C1': '#308014', 'C50c': '#49E20E', 'C66': '#C5E3BF', 'other': '#F8F8F8'}

        # Get the colours of the seqs and profiles by doing a no plot bar
        self.spb = SPBars(
            seq_count_table_path=self.sp_seq_count_path,
            profile_count_table_path=self.sp_profile_abund_and_meta_path,
            plot_type='seq_only', orientation='v', legend=False, relative_abundance=True,
            limit_genera=['C', 'D']
        )
        self.seq_color_dict = self.spb.seq_color_dict
        self.profile_color_dict = self.spb.profile_color_dict

    def _plot_cluster_leg(self, ax):
        """ Given an axis, plot a legend that contains the d and c clusterings"""
        # We will assume that this is always going to be a wide axis and we will
        # aim to plot just one row.
        # there are 2 Ds, 3 Cs, and an other.
        # so split into 6
        ax.set_xlim(0, 6)
        ax.set_ylim(0, 1)
        cluster_c_map = {'D1/D2d': '#104E8B', 'D1/D4': '#60AFFE', 'C1': '#308014', 'C50c': '#49E20E', 'C66': '#C5E3BF', 'other': '#F8F8F8'}
        rect_list = []
        for i, cluster in enumerate(['C1', 'C50c', 'C66', 'D1/D2d', 'D1/D4', 'other']):
            rect_list.append(Rectangle((i, 0.2), width=0.2, height=0.6, facecolor=cluster_c_map[cluster], edgecolor='black'))
            ax.text(x=i + 0.25, y=0.6, s=f"'{cluster}-cluster'", ha='left', va='top')
        collection = PatchCollection(rect_list, match_original=True)
        ax.add_collection(collection)
        self._rm_all_spines_and_ticks(ax)
        ax.set_xlabel('ITS2 sequence clusters')

    def _rm_all_spines_and_ticks(self, ax):
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_yticks([])
        ax.set_xticks([])

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
        seq_count_df = seq_count_df.iloc[:, index_of_first_seq:]
        self.c_sample_uid_to_post_med_absolute_dict, self.c_sample_uid_to_post_med_unique_dict, self.c_sample_uid_to_relative_genera_abund_dict = self._make_clade_specific_post_med_absolute_and_unique_dicts(clade='C', seq_count_df=seq_count_df)
        self.d_sample_uid_to_post_med_absolute_dict, self.d_sample_uid_to_post_med_unique_dict, self.d_sample_uid_to_relative_genera_abund_dict = self._make_clade_specific_post_med_absolute_and_unique_dicts(
            clade='D', seq_count_df=seq_count_df)
        return seq_count_df

    def _make_clade_specific_post_med_absolute_and_unique_dicts(self, clade, seq_count_df):
        abs_dict = {}
        unique_dict = {}
        rel_dict = {}
        for sample_ind in seq_count_df.index:
            ser = seq_count_df.loc[sample_ind][list([_ for _ in seq_count_df if (_.startswith(clade) or _.endswith(clade))])]
            ser = ser[ser != 0]
            clade_abs_abund = ser.sum()
            tot_sample_abs_abund = seq_count_df.loc[sample_ind].sum()
            abs_dict[sample_ind] = clade_abs_abund
            unique_dict[sample_ind] = len(ser)
            if clade_abs_abund == 0 or tot_sample_abs_abund == 0:
                rel_dict[sample_ind] = 0
            else:
                rel_dict[sample_ind] = ser.sum()/seq_count_df.loc[sample_ind].sum()

        return abs_dict, unique_dict, rel_dict


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
        # TODO it would make most sense to screen by relative abundance of the given clade in the sample
        # as well as the absolute and unique so we should make a dict for that too
        fig = plt.figure(figsize=(10, 10))
        # 2 down 3 across
        gs = gridspec.GridSpec(2, 3)
        axes = []
        plot_tyes = ['post_med_absolute', 'post_med_unique', 'relative']
        for i, genus in enumerate(['C', 'D']):
            for j, plot_type in enumerate(plot_tyes):
                axes.append(plt.subplot(gs[i:i+1, j:j+1]))

        # D absolute
        d_post_med_absolute_values = [self.d_sample_uid_to_post_med_absolute_dict[_] for _ in self.d_sample_uids_to_plot_non_filtered]
        # d_abs_kde = stats.gaussian_kde(d_post_med_absolute_values)
        d_abs_bins = range(0, 50000, int(50000 / 20))
        # d_abs_kde_x = np.linspace(0, 50000, 100)
        axes[0].hist(d_post_med_absolute_values, bins=d_abs_bins)
        # sec_ax = axes[0].twinx()
        # sec_ax.plot(d_abs_kde_x, d_abs_kde(d_abs_kde_x), color='black', zorder=2)
        # sec_ax.set_ylabel('Density')
        axes[0].set_ylabel('Count')
        axes[0].set_title('$\it{Durusdinium}$\npost-MED absolute sequences')
        axes[0].set_xlabel('post-MED absolute sequences')



        # D unique
        d_post_med_unique_values = [self.d_sample_uid_to_post_med_unique_dict[_] for _ in self.d_sample_uids_to_plot_non_filtered]
        # d_unique_kde = stats.gaussian_kde(d_post_med_unique_values)
        d_unique_bins = range(0, 40, 2)
        # d_unique_kde_x = np.linspace(0, 40, 100)
        axes[1].hist(d_post_med_unique_values, bins=d_unique_bins)
        # sec_ax = axes[1].twinx()
        # sec_ax.plot(d_unique_kde_x, d_unique_kde(d_unique_kde_x), color='black', zorder=2)
        # sec_ax.set_ylabel('Density')
        axes[1].set_ylabel('Count')
        axes[1].set_title('$\it{Durusdinium}$\npost-MED unique sequences')
        axes[1].set_xlabel('post-MED unique sequences')
        axes[1].ticklabel_format(useOffset=False, style='plain')

        # D relative
        d_relative_values = [self.d_sample_uid_to_relative_genera_abund_dict[_] for _ in
                                    self.d_sample_uids_to_plot_non_filtered]
        # d_rel_kde = stats.gaussian_kde(d_relative_values)
        d_rel_bins = np.arange(0, 1, 1/20)
        # d_rel_kde_x = np.linspace(0, 1, 100)
        axes[2].hist(d_relative_values, bins=d_rel_bins)
        # sec_ax = axes[2].twinx()
        # sec_ax.plot(d_rel_kde_x, d_rel_kde(d_rel_kde_x), color='black', zorder=2)
        # sec_ax.set_ylabel('Density')
        axes[2].set_ylabel('Count')
        axes[2].set_title('$\it{Durusdinium}$\nrelative abundance in sample')
        axes[2].set_xlabel('relative abundance in sample')


        # C absolute
        c_post_med_absolute_values = [self.c_sample_uid_to_post_med_absolute_dict[_] for _ in
                                        self.c_sample_uids_to_plot_non_filtered]
        # c_abs_kde = stats.gaussian_kde(c_post_med_absolute_values)
        c_abs_bins = range(0, 50000, int(50000 / 20))
        # c_abs_kde_x = np.linspace(0, 50000, 100)
        axes[3].hist(c_post_med_absolute_values, bins=c_abs_bins)
        # sec_ax = axes[3].twinx()
        # sec_ax.plot(c_abs_kde_x, c_abs_kde(c_abs_kde_x), color='black', zorder=2)
        # sec_ax.set_ylabel('Density')
        axes[3].set_ylabel('Count')
        axes[3].set_title('$\it{Cladocopium}$\npost-MED absolute sequences')
        axes[3].set_xlabel('post-MED absolute sequences')


        # C unique
        c_post_med_unique_values = [self.c_sample_uid_to_post_med_unique_dict[_] for _ in
                                  self.c_sample_uids_to_plot_non_filtered]
        # c_unique_kde = stats.gaussian_kde(c_post_med_unique_values)
        c_unique_bins = range(0, 40, 2)
        # c_unique_kde_x = np.linspace(0, 40, 100)
        axes[4].hist(c_post_med_unique_values, bins=c_unique_bins)
        # sec_ax = axes[4].twinx()
        # sec_ax.plot(c_unique_kde_x, c_unique_kde(c_unique_kde_x), color='black', zorder=2)
        # sec_ax.set_ylabel('Density')
        axes[4].set_ylabel('Count')
        axes[4].set_title('$\it{Cladocopium}$\npost-MED unique sequences')
        axes[4].set_xlabel('post-MED unique sequences')


        # C relative
        c_relative_values = [self.c_sample_uid_to_relative_genera_abund_dict[_] for _ in
                             self.c_sample_uids_to_plot_non_filtered]
        # c_rel_kde = stats.gaussian_kde(c_relative_values)
        c_rel_bins = np.arange(0, 1, 1/20)
        # c_rel_kde_x = np.linspace(0, 1, 100)
        axes[5].hist(c_relative_values, bins=c_rel_bins)
        # sec_ax = axes[5].twinx()
        # sec_ax.plot(c_rel_kde_x, c_rel_kde(c_rel_kde_x), color='black', zorder=2)
        # sec_ax.set_ylabel('Density')
        axes[5].set_ylabel('Count')
        axes[5].set_title('$\it{Cladocopium}$\nrelative abundance in sample')
        axes[5].set_xlabel('relative abundance in sample')


        plt.tight_layout()

        plt.savefig(os.path.join(self.root_dir, 'figures',
                                 f"sup_histograms_{str(datetime.now()).split('.')[0].replace('-', '').replace(' ', 'T').replace(':', '')}.svg"),
                    dpi=1200)
        plt.savefig(os.path.join(self.root_dir, 'figures',
                                 f"sup_histograms_{str(datetime.now()).split('.')[0].replace('-', '').replace(' ', 'T').replace(':', '')}.png"),
                    dpi=1200)

    def plot_supporting_hierarcical_clustering_figure(self):
        for clade in ['C', 'D']:
            self.fig = plt.figure(figsize=(10, 10))
            # 8 down 1 across
            # TODO we will need to adjust this as we refine the figure
            self.gs = gridspec.GridSpec(10, 1)
            self.axes = []
            plot_tyes = ['hierarchical', 'seq_prof', 'cluster', 'species', 'absolute', 'unique', 'relative', 'absolute_bin', 'unique_bin', 'relative_bin']

            for j, plot_type in enumerate(plot_tyes):
                self.axes.append(plt.subplot(self.gs[j: j + 1, :]))

            if clade == 'D':
                axes = [*self.axes[:4]]
                dist_output_path = self.sp_between_smp_dist_path_d
                clade_list = ['D']
                # d_clustering_dict = {(uid): (1 if rel_abund >= 0.01 else 0) for uid, rel_abund in
                #                      self.sample_uid_to_d2d_rel_abund_dict.items()}
                self._plot_for_clade(
                    axes=axes, clade_list=clade_list,
                    dist_output_path=dist_output_path,
                    sample_uids_to_plot=self.d_sample_uids_to_plot_non_filtered, cluster_dict=self.d_clustering_dict, cluster_c_map=self.d_cluster_c_map)

                foo = 'bar'
                post_med_absolute_values = [self.d_sample_uid_to_post_med_absolute_dict[_] for _ in self.d_sample_uids_to_plot_non_filtered]
                c_map = cm.get_cmap('plasma')
                self.plot_categorical_bars(ax=self.axes[4], cat_list=post_med_absolute_values, c_map=c_map)
                self.axes[4].set_xticks([])
                self.axes[4].set_yticks([])
                self.axes[4].set_ylabel('post-MED\nabsolute', fontsize='small')

                ex_absolute = [1 if _ >= 5000 else 0 for _ in post_med_absolute_values]
                filter_color_dict = {1:'white', 0:'black'}
                self.plot_categorical_bars(ax=self.axes[5], cat_list=ex_absolute, c_map=filter_color_dict)
                self.axes[5].set_xticks([])
                self.axes[5].set_yticks([])
                self.axes[5].set_ylabel('post-MED\nabsolute', fontsize='small')

                post_med_unique_values = [self.d_sample_uid_to_post_med_unique_dict[_] for _ in self.d_sample_uids_to_plot_non_filtered]
                c_map = cm.get_cmap('plasma')
                self.plot_categorical_bars(ax=self.axes[6], cat_list=post_med_unique_values, c_map=c_map)
                self.axes[6].set_xticks([])
                self.axes[6].set_yticks([])
                self.axes[6].set_ylabel('post-MED\nunique', fontsize='small')

                ex_unique = [1 if _ >= 10 else 0 for _ in post_med_unique_values]
                self.plot_categorical_bars(ax=self.axes[7], cat_list=ex_unique, c_map=filter_color_dict)
                self.axes[7].set_xticks([])
                self.axes[7].set_yticks([])
                self.axes[7].set_ylabel('post-MED\nunique', fontsize='small')

                post_med_relative_values = [self.d_sample_uid_to_relative_genera_abund_dict[_] for _ in
                                          self.d_sample_uids_to_plot_non_filtered]
                c_map = cm.get_cmap('plasma')
                self.plot_categorical_bars(ax=self.axes[8], cat_list=post_med_relative_values, c_map=c_map)
                self.axes[8].set_xticks([])
                self.axes[8].set_yticks([])
                self.axes[8].set_ylabel('post-MED\nrelative abund.', fontsize='small')
                ex_relative = [1 if _ >= 0.25 else 0 for _ in post_med_relative_values]
                self.plot_categorical_bars(ax=self.axes[9], cat_list=ex_relative, c_map=filter_color_dict)
                self.axes[9].set_xticks([])
                self.axes[9].set_yticks([])
                self.axes[9].set_ylabel('post-MED\nrelative abund.', fontsize='small')

            elif clade == 'C':
                axes = [*self.axes[:4]]
                dist_output_path = self.sp_between_smp_dist_path_c
                clade_list = ['C']
                cluster_to_number_map = {'C1': 0, 'C50c': 0.25, 'C66': 0.5, 'other': 1}
                # c_clustering_dict = {uid: cluster_to_number_map[self.c_clustering_dict[uid]] for uid in
                #                      self.c_sample_uids_to_plot_non_filtered}
                self._plot_for_clade(
                    axes=axes, clade_list=clade_list,
                    dist_output_path=dist_output_path,
                    sample_uids_to_plot=self.c_sample_uids_to_plot_non_filtered,
                    cluster_dict=self.c_clustering_dict, cluster_c_map=self.c_cluster_c_map)

                foo = 'bar'
                c_map = cm.get_cmap('plasma')
                filter_color_dict = {1: 'white', 0: 'black'}
                post_med_absolute_values = [self.c_sample_uid_to_post_med_absolute_dict[_] for _ in
                                            self.c_sample_uids_to_plot_non_filtered]
                self.plot_categorical_bars(ax=self.axes[4], cat_list=post_med_absolute_values, c_map=c_map)
                self.axes[4].set_xticks([])
                self.axes[4].set_yticks([])
                self.axes[4].set_ylabel('post-MED\nabsolute', fontsize='small')
                ex_absolute = [1 if _ >= 5000 else 0 for _ in post_med_absolute_values]
                self.plot_categorical_bars(ax=self.axes[5], cat_list=ex_absolute, c_map=filter_color_dict)
                self.axes[5].set_xticks([])
                self.axes[5].set_yticks([])
                self.axes[5].set_ylabel('post-MED\nabsolute', fontsize='small')

                post_med_unique_values = [self.c_sample_uid_to_post_med_unique_dict[_] for _ in
                                          self.c_sample_uids_to_plot_non_filtered]
                self.plot_categorical_bars(ax=self.axes[6], cat_list=post_med_unique_values, c_map=c_map)
                self.axes[6].set_xticks([])
                self.axes[6].set_yticks([])
                self.axes[6].set_ylabel('post-MED\nunique', fontsize='small')
                ex_unique = [1 if _ >= 10 else 0 for _ in post_med_unique_values]
                self.plot_categorical_bars(ax=self.axes[7], cat_list=ex_unique, c_map=filter_color_dict)
                self.axes[7].set_xticks([])
                self.axes[7].set_yticks([])
                self.axes[7].set_ylabel('post-MED\nunique', fontsize='small')

                post_med_relative_values = [self.c_sample_uid_to_relative_genera_abund_dict[_] for _ in
                                            self.c_sample_uids_to_plot_non_filtered]
                self.plot_categorical_bars(ax=self.axes[8], cat_list=post_med_relative_values, c_map=c_map)
                self.axes[8].set_xticks([])
                self.axes[8].set_yticks([])
                self.axes[8].set_ylabel('post-MED\nrelative abund.', fontsize='small')
                ex_relative = [1 if _ >= 0.25 else 0 for _ in post_med_relative_values]
                self.plot_categorical_bars(ax=self.axes[9], cat_list=ex_relative, c_map=filter_color_dict)
                self.axes[9].set_xticks([])
                self.axes[9].set_yticks([])
                self.axes[9].set_ylabel('post-MED\nrelative abund.', fontsize='small')

            plt.tight_layout()
            print('saving svg')
            plt.savefig(os.path.join(self.root_dir, 'figures',
                                     f"supporting_hierarchical_clade_{clade}_{str(datetime.now()).split('.')[0].replace('-', '').replace(' ', 'T').replace(':', '')}.svg"),
                        dpi=1200)
            print('saving png')
            plt.savefig(os.path.join(self.root_dir, 'figures',
                                     f"supporting_hierarchical_clade_{clade}_{str(datetime.now()).split('.')[0].replace('-', '').replace(' ', 'T').replace(':', '')}.png"),
                        dpi=1200)
            foo = 'bar'

    def plot_main_hierarchical_clutering_figure(self):
        # Get the list of samples that need plotting
        # self.sp_between_smp_dist_path_c
        # self.sp_between_smp_dist_path_d
        self.fig = plt.figure(figsize=(10, 10))
        # 6 down 4 across
        # TODO we will need to adjust this as we refine the figure
        self.gs = gridspec.GridSpec(8, 3)
        self.axes = []
        plot_tyes = ['hierarchical', 'seq_prof', 'cluster', 'species']
        for i, genus in enumerate(['Durusdinium', 'Cladocopium']):
            for j, plot_type in enumerate(plot_tyes):
                if i == 0:
                    self.axes.append(plt.subplot(self.gs[((i * len(plot_tyes)) + j):((i * len(plot_tyes)) + j) + 1, :]))
                else:
                    self.axes.append(plt.subplot(self.gs[((i * len(plot_tyes)) + j):((i * len(plot_tyes)) + j) + 1, :2]))
        # Legend axis
        self.leg_ax = plt.subplot(self.gs[4:, 2:])
        self.spb.plot_only_legend(seq_leg_ax=self.leg_ax)
        self.leg_ax.set_title('20 most abundant\nITS2 sequences')
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
        # Plot D
        d_samples_to_plot = [uid for uid in self.d_sample_uids_to_plot_non_filtered if (
                (self.d_sample_uid_to_post_med_absolute_dict[uid] >= 5000) and
                (self.d_sample_uid_to_post_med_unique_dict[uid] >= 10)
        )]
        print(f"{len(self.d_sample_uids_to_plot_non_filtered)-len(d_samples_to_plot)} D samples removed by filtering")
        axes = [*self.axes[:4]]
        dist_output_path = self.sp_between_smp_dist_path_d
        clade_list = ['D']
        self._plot_for_clade(
            axes=axes, clade_list=clade_list,
            dist_output_path=dist_output_path,
            sample_uids_to_plot=d_samples_to_plot, cluster_dict=self.d_clustering_dict, cluster_c_map=self.d_cluster_c_map)


        foo  ='bar'
        c_samples_to_plot = [uid for uid in self.c_sample_uids_to_plot_non_filtered if (
                (self.c_sample_uid_to_post_med_absolute_dict[uid] >= 5000) and
                (self.c_sample_uid_to_post_med_unique_dict[uid] >= 10)
        )]
        print(f"{len(self.c_sample_uids_to_plot_non_filtered) - len(c_samples_to_plot)} C samples removed by filtering")
        axes = [*self.axes[4:]]
        dist_output_path = self.sp_between_smp_dist_path_c
        clade_list = ['C']

        # c_clustering_dict = {uid: self.c_clustering_dict[uid] for uid in
        #                      self.c_sample_uids_to_plot_non_filtered}
        self._plot_for_clade(
            axes=axes, clade_list=clade_list,
            dist_output_path=dist_output_path, sample_uids_to_plot=c_samples_to_plot,
            cluster_dict=self.c_clustering_dict, cluster_c_map=self.c_cluster_c_map
        )

        foo = 'bar'
        plt.tight_layout()
        print('saving svg')
        plt.savefig(os.path.join(self.root_dir, 'figures',
                                 f"main_hierarchical_clustering_{str(datetime.now()).split('.')[0].replace('-', '').replace(' ', 'T').replace(':', '')}.svg"),
                    dpi=1200)
        print('saving png')
        plt.savefig(os.path.join(self.root_dir, 'figures',
                                 f"main_hierarchical_clustering_{str(datetime.now()).split('.')[0].replace('-', '').replace(' ', 'T').replace(':', '')}.png"),
                    dpi=1200)

    def switch_off_spines(self, ax):
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

    def _plot_for_clade(self, axes, clade_list, dist_output_path, sample_uids_to_plot, cluster_dict, cluster_c_map):

        sph_plot = SPHierarchical(dist_output_path=dist_output_path, ax=axes[0],
                                  sample_uids_included=sample_uids_to_plot)
        sph_plot.plot()
        # Thin out the lines
        axes[0].collections[0].set_linewidth(0.5)
        self.switch_off_spines(ax=axes[0])
        dendrogram_sample_uid_order = sph_plot.dendrogram['ivl']
        # Now plot up the sequences and profiles
        spb = SPBars(
            seq_count_table_path=self.sp_seq_count_path,
            profile_count_table_path=self.sp_profile_abund_and_meta_path,
            plot_type='seq_only', orientation='h', legend=False, relative_abundance=True,
            sample_uids_included=dendrogram_sample_uid_order, bar_ax=axes[1], limit_genera=clade_list,
            seq_profile_scalar=(1.0, 0.3), seq_color_dict=self.seq_color_dict,
            profile_color_dict=self.profile_color_dict
        )
        spb.plot()
        axes[1].set_xticks([])
        axes[1].set_yticks([])
        cluster_vals = [cluster_dict[uid] for uid in dendrogram_sample_uid_order]
        self.plot_categorical_bars(ax=axes[2], cat_list=cluster_vals, c_map=cluster_c_map)
        species_vals = [self.sp_datasheet_df.at[self.sp_sample_uid_to_sample_name_dict[_], 'host_species'] for _ in dendrogram_sample_uid_order]

        self.plot_categorical_bars(ax=axes[3], cat_list=species_vals, c_map=self.species_c_map)

        if 'C' in clade_list:
            axes[0].set_title('Cladocopium', fontsize='medium', style='italic')
        else:
            axes[0].set_title('Durusdinium', fontsize='medium', style='italic')
        for ax in axes:
            ax.set_xticks([])
            ax.set_yticks([])

        axes[0].set_ylabel('BrayCurtis\ndissimilarity', fontsize='small')
        axes[1].set_ylabel('ITS2 sequence\ndiversity', fontsize='small')
        axes[2].set_ylabel('assigned\ncluster', fontsize='small')
        axes[3].set_ylabel('host\nspecies', fontsize='small')

    def plot_categorical_bars(self, ax, cat_list, c_map):
        # TODO accept an acutal colour map and use this for the rectangles
        left = 0
        rect_list = []
        if type(c_map) == mpl.colors.ListedColormap:
            # normalise the cat_list
            cat_list = np.array(cat_list).reshape(-1, 1)
            min_max_scaler = MinMaxScaler()
            cat_list_norm = min_max_scaler.fit_transform(cat_list)
            # Working with a custom colour dict
            for cat_item in cat_list_norm:
                rect_list.append(Rectangle(
                    (left, 0),
                    1,
                    1, color=c_map(cat_item[0])))
                left += 1
        else:
            # Working with a colour dict
            for cat_item in cat_list:
                rect_list.append(Rectangle(
                    (left, 0),
                    1,
                    1, color=c_map[cat_item]))
                left += 1
        patches_collection = PatchCollection(rect_list, match_original=True)
        ax.add_collection(patches_collection)
        ax.set_xlim(0, len(cat_list))
        ax.set_ylim(0, 1)

class ClusteredZooxs(SchuppFigures):
    def __init__(self, color_by_cluster=True):
        """
        The base for the main figure showing the zooxs results.
        Four columns, one per species, and two rows, adults and recruits.
        The recruits will be a matrix of time and temperature

        If color_by_cluster is set to False then we will plot the full sequence colors
        """
        super().__init__()
        self.color_by_cluster = color_by_cluster
        self.species_full = [
            'Acropora digitifera', 'Acropora hyacinthus', 'Pocillopora damicornis', 'Leptastrea purpurea'
        ]
        self.species_short = ['ad', 'ah', 'pd', 'lp']
        self.species_short_to_full_dict = {
            'ad': 'Acropora digitifera', 'ah': 'Acropora hyacinthus',
            'pd': 'Pocillopora damicornis', 'lp': 'Leptastrea purpurea'
        }
        self.fig = plt.figure(figsize=(15, 5))
        self.gs = gridspec.GridSpec(5, 4)
        self.species_axes_dict = {}
        for i, species in enumerate(self.species_short):
            temp_list = []
            # First the adult axis
            temp_list.append(plt.subplot(self.gs[:1, i]))
            # Then the recruit matrix
            inner_grid_spec = self.gs[1:4, i].subgridspec(3, 5)
            outer_temp_list = []
            for k in range(3):
                inner_temp_list = []
                for l in range(5):
                    inner_temp_list.append(plt.subplot(inner_grid_spec[k, l]))
                outer_temp_list.append(inner_temp_list)
            temp_list.append(outer_temp_list)
            self.species_axes_dict[species] = temp_list
        self.leg_ax = plt.subplot(self.gs[4:,:])
        if self.color_by_cluster:
            self._plot_cluster_leg(ax=self.leg_ax)
        else:
            self.spb.plot_only_legend(seq_leg_ax=self.leg_ax)
            self.leg_ax.set_xlabel('Top 20 most abundant ITS2 sequences')
        foo = 'bar'

    def plot(self):
        for sp in self.species_short:
            self._plot_adult_zooxs(sp=sp, ax=self.species_axes_dict[sp][0])
            self._plot_recruit_zooxs(sp=sp, ax_array=self.species_axes_dict[sp][1])
        plt.tight_layout()
        if self.color_by_cluster:
            print('saving svg')
            plt.savefig(os.path.join(self.root_dir, 'figures',
                                     f"zooxs_base_figure_cluster_{str(datetime.now()).split('.')[0].replace('-', '').replace(' ', 'T').replace(':', '')}.svg"),
                        dpi=1200)
            print('saving png')
            plt.savefig(os.path.join(self.root_dir, 'figures',
                                     f"zooxs_base_figure_cluster_{str(datetime.now()).split('.')[0].replace('-', '').replace(' ', 'T').replace(':', '')}.png"),
                        dpi=1200)
        else:
            print('saving svg')
            plt.savefig(os.path.join(self.root_dir, 'figures',
                                     f"zooxs_base_figure_no_cluster_{str(datetime.now()).split('.')[0].replace('-', '').replace(' ', 'T').replace(':', '')}.svg"),
                        dpi=1200)
            print('saving png')
            plt.savefig(os.path.join(self.root_dir, 'figures',
                                     f"zooxs_base_figure_no_cluster_{str(datetime.now()).split('.')[0].replace('-', '').replace(' ', 'T').replace(':', '')}.png"),
                        dpi=1200)

    def _plot_temp_time_recruit_zooxs_matrix(self, ax_array, working_datasheet_df, sp):
        # Now we want to plot up the rectangles on a per temperature/time point combinations basis.
        # We should be able to use the rectangle code that we made for the adults
        # To input into that code we simply need a list of sample UIDs and an axis
        for k, temp in enumerate([29, 30, 31]):
            for l, age in enumerate([1, 3, 6, 9, 12]):
                ax = ax_array[k][l]
                # If left hand plot, set temp as y axis label
                # If bottom plot, set time as x axis label
                if k ==0 and l == 2:
                    ax.set_title('Recruit', fontsize='small')
                if l == 0 and k == 1 and sp == 'ad':
                    ax.set_ylabel(f'Temperature\n{temp} {DEGREE_SIGN}C', fontsize='small')
                elif l == 0 and sp == 'ad':
                    ax.set_ylabel(f'{temp} {DEGREE_SIGN}C', fontsize='small')
                if k == 2 and l == 2:
                    # We will try to put month in the middle
                    ax.set_xlabel(f'{age}\nmonths')
                elif k==2:
                    ax.set_xlabel(f'{age}')
                sample_uids = [
                    self.sp_sample_name_to_sample_uid_dict[sample_name] for
                    sample_name in
                    working_datasheet_df[
                        (working_datasheet_df['temp'] == temp) &
                        (working_datasheet_df['age'] == age)
                        ].index
                ]
                if sample_uids:
                    ax.set_xticks([])
                    ax.set_yticks([])
                    # Not all species have the zooxs data for the complete time/temp matrix
                    self._plot_seq_rectangles_adult_zooxs(ax=ax, sample_uids=sample_uids)
                else:
                    ax.set_xticks([])
                    ax.set_yticks([])
                    ax.set_ylim(0, 1)
                    ax.set_xlim(0, 1)
                    ax.text(0.5, 0.5, 'no\nsamples', ha='center', va='center', fontsize='xx-small')


    def _plot_recruit_zooxs(self, sp, ax_array):
        sample_uids, sample_names_index = self._get_sample_uid_recruit_zooxs(sp)
        working_datasheet_df = self.sp_datasheet_df.loc[sample_names_index]
        # convert the ages to numeric by getting rid of the 'month' or 'months'
        working_datasheet_df['age'] = [int(v.split(' ')[0]) for k, v in working_datasheet_df['age'].iteritems()]
        self._plot_temp_time_recruit_zooxs_matrix(ax_array, working_datasheet_df, sp)

    def _plot_adult_zooxs(self, sp, ax):
        sample_uids = self._get_sample_uid_adult_zooxs(sp)
        # Now we can plot up the seqs and profiles.
        self._plot_seq_rectangles_adult_zooxs(ax, sample_uids)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f'$\it{self.species_short_to_full_dict[sp]}$\nAdult', fontsize='small')

    def _plot_seq_rectangles_adult_zooxs(self, ax, sample_uids):
        """"""
        x_index_for_plot = 0
        patches_list = []
        if self.color_by_cluster:
            # Then we want to plot color by cluster
            # To do this, for each sample, look to see if there is C and or D in each sample
            # Then for each of C and D, see if we can look them up in the respective
            # cluster dict. If we can, then plot that colour else plot other color
            for sample_uid in sample_uids:
                bottom = 0
                for g in ['D', 'C']:
                    g_seqs = self.sp_seq_rel_abund_df.loc[
                        sample_uid,
                        [_ for _ in list(self.sp_seq_rel_abund_df) if (_.startswith(g) or _.endswith(g))]
                    ]
                    g_seqs_tot = g_seqs.sum()
                    if g_seqs_tot > 0:
                        # Then there are genus seqs in this sample
                        if g == 'C':
                            if sample_uid in self.c_clustering_dict:
                                # Then we can assign a cluster
                                color = self.c_cluster_c_map[self.c_clustering_dict[sample_uid]]
                            else:
                                # Assign the 'other' color for the genus
                                color = self.c_cluster_c_map['other']
                        else:
                            if sample_uid in self.d_clustering_dict:
                                # Then we can assign a cluster
                                color = self.d_cluster_c_map[self.d_clustering_dict[sample_uid]]
                            else:
                                # Assign the 'other' color for the genus
                                color = self.d_cluster_c_map['other']
                        patches_list.append(Rectangle((x_index_for_plot,bottom), width=1, height=g_seqs_tot, color=color))
                    bottom += g_seqs_tot
                if x_index_for_plot != 0:
                    ax.axvline(x_index_for_plot, zorder=2, color='black', linewidth=0.5)
                x_index_for_plot += 1

                # we will also want to add vertical lines to denote the samples
            patches_collection = PatchCollection(patches_list, match_original=True, zorder=1)
            ax.add_collection(patches_collection)
            ax.set_xlim(0, len(sample_uids))
            ax.set_ylim(0, 1)
        else:
            spb = self.spb = SPBars(
            seq_count_table_path=self.sp_seq_count_path,
            plot_type='seq_only', orientation='h', legend=False, relative_abundance=True,
            limit_genera=['C', 'D'], sample_uids_included=sample_uids, bar_ax=ax, seq_color_dict=self.seq_color_dict
            )
            spb.plot()
            # for sample_uid in sample_uids:
            #     bottom = 0
            #     non_zero_seq_abundances = self.sp_seq_rel_abund_df.loc[sample_uid][
            #         self.sp_seq_rel_abund_df.loc[sample_uid] > 0]
            #     for seq_uid, rel_abund in non_zero_seq_abundances.iteritems():
            #         patches_list.append(Rectangle(
            #             (x_index_for_plot - 0.5, bottom),
            #             1,
            #             rel_abund, color=self.seq_color_dict[seq_uid]))
            #         bottom += rel_abund
            #     x_index_for_plot += 1


    def _get_sample_uid_recruit_zooxs(self, sp):
        if sp == 'ad':
            # 1, 3, 6, 9, 12
            # 29, 30, 31
            sample_names_index = self.sp_datasheet_df[
                    (self.sp_datasheet_df['host_species'] == 'digitifera') &
                    (self.sp_datasheet_df['age'].str.contains("month")) &
                    (~self.sp_datasheet_df['age'].str.contains("through")) &
                    (~self.sp_datasheet_df['age'].str.contains("delayed"))
                    ].index
            sample_uids = [
                self.sp_sample_name_to_sample_uid_dict[sample_name] for sample_name in sample_names_index
            ]
        elif sp == 'ah':
            # 3, 12
            # 29, 30, 31
            sample_names_index = self.sp_datasheet_df[
                (self.sp_datasheet_df['host_species'] == 'surculosa') &
                (self.sp_datasheet_df['age'].str.contains("month")) &
                (~self.sp_datasheet_df['age'].str.contains("adult"))
                ].index
            sample_uids = [
                self.sp_sample_name_to_sample_uid_dict[sample_name] for sample_name in sample_names_index
            ]
        elif sp == 'pd':
            # 12
            # 29, 30, 31
            sample_names_index = self.sp_datasheet_df[
                (self.sp_datasheet_df['host_species'] == 'damicornis') &
                (self.sp_datasheet_df['age'].str.contains("month")) &
                (~self.sp_datasheet_df['age'].str.contains("adult"))
                ].index
            sample_uids = [
                self.sp_sample_name_to_sample_uid_dict[sample_name] for sample_name in sample_names_index
            ]
        elif sp == 'lp':
            # 12
            # 29, 30, 31
            sample_names_index = self.sp_datasheet_df[
                (self.sp_datasheet_df['host_species'] == 'purpurea') &
                (self.sp_datasheet_df['age'].str.contains("month")) &
                (~self.sp_datasheet_df['age'].str.contains("adult"))
                ].index
            sample_uids = [
                self.sp_sample_name_to_sample_uid_dict[sample_name] for sample_name in sample_names_index
            ]
        else:
            raise RuntimeError(f'unexpected species {sp}')
        return sample_uids, sample_names_index

    def _get_sample_uid_adult_zooxs(self, sp):
        if sp == 'ad':
            sample_uids = [
                self.sp_sample_name_to_sample_uid_dict[sample_name] for sample_name in self.sp_datasheet_df[
                    (self.sp_datasheet_df['age'] == 'adult') &
                    (self.sp_datasheet_df['host_species'] == 'digitifera')
                    ].index
            ]
        elif sp == 'ah':
            sample_uids = [
                self.sp_sample_name_to_sample_uid_dict[sample_name] for sample_name in self.sp_datasheet_df[
                    (self.sp_datasheet_df['age'] == 'adult') &
                    (self.sp_datasheet_df['host_species'] == 'surculosa')
                    ].index
            ]
        elif sp == 'pd':
            sample_uids = [
                self.sp_sample_name_to_sample_uid_dict[sample_name] for sample_name in self.sp_datasheet_df[
                    (self.sp_datasheet_df['age'] == 'adult') &
                    (self.sp_datasheet_df['host_species'] == 'damicornis')
                    ].index
            ]
        elif sp == 'lp':
            # NB there are multiple purpurea adult samples from different location
            # I'm not sure which onces I should be working with
            # For the time being I'll proceed with thos samples that were part of the same sequencing
            # plate as the other adults
            sample_uids = [
                self.sp_sample_name_to_sample_uid_dict[sample_name] for sample_name in self.sp_datasheet_df[
                    (self.sp_datasheet_df['age'] == 'adult') &
                    (self.sp_datasheet_df['host_species'] == 'purpurea')
                    ].index if
                "P4" in sample_name
            ]
        else:
            raise RuntimeError(f'unexpected species {sp}')
        return sample_uids

h = HierarchicalPlot()
h.plot_main_hierarchical_clutering_figure()
# h.plot_supporting_hierarcical_clustering_figure()
# h.plot_supporting_histograms()
# ClusteredZooxs(color_by_cluster=False).plot()
# ClusteredZooxs(color_by_cluster=True).plot()

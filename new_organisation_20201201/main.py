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
import itertools
import sys
import random
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
DEGREE_SIGN = u'\N{DEGREE SIGN}'


class schupp_figures:
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
        self.sp_profile_abund_path = os.path.join(
            self.sp_data_path, 'its2_type_profiles',
            '125_20201018_DBV_20201020T020625.profiles.absolute.abund_only.txt'
        )
        self.sp_profile_meta_path = os.path.join(
            self.sp_data_path, 'its2_type_profiles',
            '125_20201018_DBV_20201020T020625.profiles.meta_only.txt'
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
        self.data_types = ['adult_survival', 'adult_zooxs', 'recruit_survival', 'recruit_size', 'recruit_fv_fm', 'recruit_zooxs']

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
        self.sp_profile_rel_abund_df = self.sp_profile_abs_abund_df.div(
            self.sp_profile_abs_abund_df.sum(axis=1), axis=0
        )

        # Reorder the abundance dataframes in order of most abundant sequences/profiles
        sorted_seq_index = self.sp_seq_rel_abund_df.sum(axis=0).sort_values(ascending=False).index
        sorted_profile_index = self.sp_profile_rel_abund_df.sum(axis=0).sort_values(ascending=False).index
        self.sp_seq_rel_abund_df.reindex(sorted_seq_index, axis=1)
        self.sp_seq_abs_abund_df.reindex(sorted_seq_index, axis=1)
        self.sp_profile_rel_abund_df.reindex(sorted_profile_index, axis=1)
        self.sp_profile_abs_abund_df.reindex(sorted_profile_index, axis=1)

        # Figure
        # 10 wide, 6 deep
        self.fig = plt.figure(figsize=(10, 10))
        self.gs = gridspec.GridSpec(6, 4)
        self.axes = [[] for sp in self.species_short]
        # We will have the various axes as coordinates in columns and rows
        # For all of the data types except recruit zooxs, there will be one axis per speceis, data type combination
        # For the recruit_zooxs this gets more complicated as we esentially have a matrix of time and temperature
        # To represent this, for the recruit zooxs we will have a matrix of axes
        # We will house these axes array (one for each species) in a dict that is key of the speices short
        # and value of a 2d list
        self.recruit_zooxs_axes_dict = {}
        for i, sp in enumerate(self.species_short):
            for j, data_type in enumerate(self.data_types):
                # create a plot
                if data_type == 'recruit_zooxs':
                    # We want to split up the recruit_zooxs gs to make a matrix of plots
                    inner_grid_spec = self.gs[j, i].subgridspec(3, 5)
                    outer_temp_list = []
                    for k in range(3):
                        inner_temp_list = []
                        for l in range(5):
                            inner_temp_list.append(plt.subplot(inner_grid_spec[k, l]))
                        outer_temp_list.append(inner_temp_list)
                    self.recruit_zooxs_axes_dict[sp] = outer_temp_list
                else:
                    self.axes[i].append(plt.subplot(self.gs[j, i]))

        self.temp_plot_marker_dict = {29: "v", 30: "o", 31: "^", '29': "v", '30': "o", '31': "^"}
        self.temp_plot_colour_dict = {29: "#a6a6a6", 30: "#898989", 31: "#0d0d0d", '29': "#a6a6a6", '30': "#898989", '31': "#0d0d0d"}

        # Colour generators for seq and profile plotting
        self.grey_iterator = itertools.cycle(['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F'])
        self.colour_hash_iterator = iter(self._get_colour_list())
        self.pre_def_seq_colour_dict = self._get_pre_def_colour_dict()
        self.colour_palette_pas_gen = ('#%02x%02x%02x' % rgb_tup for rgb_tup in
                                  self._create_colour_list(mix_col=(255, 255, 255), sq_dist_cutoff=5000, num_cols=8,
                                                           time_out_iterations=10000))
        self.seq_color_dict = self._make_seq_colour_dict()
        self.prof_color_dict = self._make_profile_color_dict()

    def _make_profile_color_dict(self):
        prof_color_dict = {}
        for prof_uid in list(self.sp_profile_rel_abund_df):
            try:
                prof_color_dict[prof_uid] = next(self.colour_palette_pas_gen)
            except StopIteration:
                prof_color_dict[prof_uid] = next(self.grey_iterator)
        return prof_color_dict

    def _make_seq_colour_dict(self):
        seq_color_dict = {}
        for seq_name in list(self.sp_seq_rel_abund_df):
            if seq_name in self.pre_def_seq_colour_dict:
                seq_color_dict[seq_name] = self.pre_def_seq_colour_dict[seq_name]
            else:
                try:
                    seq_color_dict[seq_name] = next(self.colour_hash_iterator)
                except StopIteration:
                    seq_color_dict[seq_name] = next(self.grey_iterator)
        return seq_color_dict

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
        self.sp_profile_name_to_profile_uid_dict = dict(zip(
            profile_meta_df["ITS2 type profile UID"],
            profile_meta_df['ITS2 type profile']
        ))
        self.sp_profile_uid_to_profile_name_dict = dict(zip(
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
        self.sp_sample_uid_to_sample_name_dict = dict(zip(seq_count_df.sample_uid, seq_count_df.sample_name))
        self.sp_sample_name_to_sample_uid_dict = dict(zip(seq_count_df.sample_name, seq_count_df.sample_uid))
        seq_count_df.set_index(keys='sample_uid', drop=True, inplace=True)
        index_of_first_seq = list(seq_count_df).index("A3")
        seq_count_df = seq_count_df.iloc[:-1, index_of_first_seq:]
        return seq_count_df

    def _make_sp_datasheet_df(self):
        sp_ds_df = pd.read_excel(io=self.sp_datasheet_path, skiprows=[0])
        collection_cols = [lab for lab in list(sp_ds_df) if ("collection" in lab) or ("Unnamed" in lab)]
        sp_ds_df.drop(columns=collection_cols, inplace=True)
        sp_ds_df.set_index("sample_name", drop=True, inplace=True)
        return sp_ds_df

    def plot_fig_1(self):
        """
        The main plotting method for plotting up the figure that has the four species columns with a plot
        for each of the data types.
        """
        for i, sp in enumerate(self.species_short):
            for j, data_type in enumerate(self.data_types):
                if data_type == 'adult_survival':
                    self._plot_adult_survival(i, j, sp, data_type)
                if data_type == 'adult_zooxs':
                    self._plot_adult_zooxs(i, j, sp)
                elif data_type == 'recruit_survival':
                    self._plot_recruit_survival(data_type, i, j, sp)
                elif data_type == 'recruit_size':
                    self._plot_recruit_size(data_type, i, j, sp)
                elif data_type == 'recruit_fv_fm':
                    self._plot_recruit_fv_fm(data_type, i, j, sp)
                elif data_type == 'recruit_zooxs':
                    self._plot_recruit_zooxs(sp)
        foo = 'bar'

    def _plot_adult_zooxs(self, i, j, sp):
        ax = self.axes[i][j]
        sample_uids = self._get_sample_uid_adult_zooxs(sp)
        # Now we can plot up the seqs and profiles.
        self._plot_seq_rectangles_adult_zooxs(ax, sample_uids)

    def _plot_recruit_zooxs(self, sp):
        ax_array = self.recruit_zooxs_axes_dict[sp]
        sample_uids, sample_names_index = self._get_sample_uid_recruit_zooxs(sp)
        # convert the ages to numeric by getting rid of the 'month' or 'months'
        working_datasheet_df = self.sp_datasheet_df.loc[sample_names_index]
        working_datasheet_df['age'] = [int(v.split(' ')[0]) for k, v in working_datasheet_df['age'].iteritems()]
        self._plot_temp_time_recruit_zooxs_matrix(ax_array, working_datasheet_df)

    def _plot_temp_time_recruit_zooxs_matrix(self, ax_array, working_datasheet_df):
        # Now we want to plot up the rectanges on a per temperature/time point combinations basis.
        # We should be able to use the rectangle code that we made for the adults
        # To input into that code we simply need a list of sample UIDs and an axis
        for k, temp in enumerate([29, 30, 31]):
            for l, age in enumerate([1, 3, 6, 9, 12]):
                ax = ax_array[k][l]
                sample_uids = [
                    self.sp_sample_name_to_sample_uid_dict[sample_name] for
                    sample_name in
                    working_datasheet_df[
                        (working_datasheet_df['temp'] == temp) &
                        (working_datasheet_df['age'] == age)
                        ].index
                ]
                if sample_uids:
                    # Not all species have the zooxs data for the complete time/temp matrix
                    self._plot_seq_rectangles_adult_zooxs(ax=ax, sample_uids=sample_uids)

    def _plot_seq_rectangles_adult_zooxs(self, ax, sample_uids):
        x_index_for_plot = 0
        patches_list = []
        colour_list = []
        for sample_uid in sample_uids:
            bottom = 0
            non_zero_seq_abundances = self.sp_seq_rel_abund_df.loc[sample_uid][
                self.sp_seq_rel_abund_df.loc[sample_uid] > 0]
            for seq_uid, rel_abund in non_zero_seq_abundances.iteritems():
                patches_list.append(Rectangle(
                    (x_index_for_plot - 0.5, bottom),
                    1,
                    rel_abund, color=self.seq_color_dict[seq_uid]))
                bottom += rel_abund
                colour_list.append(self.seq_color_dict[seq_uid])
            x_index_for_plot += 1
        listed_colour_map = ListedColormap(colour_list)
        patches_collection = PatchCollection(patches_list, cmap=listed_colour_map)
        patches_collection.set_array(np.arange(len(patches_list)))
        ax.add_collection(patches_collection)
        ax.autoscale_view()
        self.fig.canvas.draw()

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
            # TODO there are multiple purpurea adult samples from different location
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

    def _plot_recruit_fv_fm(self, data_type, i, j, sp):
        ax = self.axes[i][j]
        param_to_plot = "fv_fm"
        working_df = self.fv_fm_size_df[
            (self.fv_fm_size_df['species'] == sp) & (self.fv_fm_size_df[param_to_plot].notnull())]
        if working_df.index.to_list():
            self._plot_a_set_of_line_data(ax, data_type, sp, working_df, param_to_plot)

    def _plot_recruit_size(self, data_type, i, j, sp):
        ax = self.axes[i][j]
        param_to_plot = "cyl_vol"
        working_df = self.fv_fm_size_df[
            (self.fv_fm_size_df['species'] == sp) & (self.fv_fm_size_df[param_to_plot].notnull())]
        # Convert to cm3
        working_df['cyl_vol'] = working_df['cyl_vol'] / 1000
        if working_df.index.to_list():
            self._plot_a_set_of_line_data(ax, data_type, sp, working_df, param_to_plot)

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
        ]
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
        if working_df.index.to_list():
            # NB for the recruit survival points, we don't want there to be a line connecting the before
            # and after travel points.
            # TODO we may want to add in a vertical line or something similar to denote the travel time for the samples
            self._plot_a_set_of_line_data(ax, data_type, sp, working_df, param_to_plot, time_first_after_val)

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
        ax.legend(loc='lower left', fontsize='xx-small')
        ax.set_xlabel(time_unit, fontsize='xx-small')
        if sp == 'ad' and "survival" in data_type:
            # only set one per row
            ax.set_ylabel('survial %', fontsize='xx-small')
        elif sp == 'ad' and "size" in data_type:
            ax.set_ylabel('cyl. vol. ml', fontsize='xx-small')
        elif sp == 'ad' and data_type == "recruit_fv_fm":
            ax.set_ylabel('Fv/Fm', fontsize='xx-small')
        if data_type == 'adult_survival':
            ax.set_title(self.species_short_to_full_dict[sp], fontsize='small')
        plt.tight_layout()

    def _calc_mean_sem_plot_line(self, ax, param_to_plot, ser, temp):
        means = [ser[ser['time_value'] == time_val][param_to_plot].mean() for time_val in
                 ser['time_value'].unique()]
        sem = [ser[ser['time_value'] == time_val][param_to_plot].sem() for time_val in
               ser['time_value'].unique()]
        ax.errorbar(
            x=ser['time_value'].unique(), y=means, yerr=sem, marker='o',
            linestyle='--', linewidth=1, ecolor=self.temp_plot_colour_dict[temp],
            elinewidth=1, color=self.temp_plot_colour_dict[temp], markersize=2,
            label=f'{temp}{DEGREE_SIGN}C'
        )

    def _calc_mean_sem_plot_line_no_label(self, ax, param_to_plot, ser, temp):
        means = [ser[ser['time_value'] == time_val][param_to_plot].mean() for time_val in
                 ser['time_value'].unique()]
        sem = [ser[ser['time_value'] == time_val][param_to_plot].sem() for time_val in
               ser['time_value'].unique()]
        ax.errorbar(
            x=ser['time_value'].unique(), y=means, yerr=sem, marker='o',
            linestyle='--', linewidth=1, ecolor=self.temp_plot_colour_dict[temp],
            elinewidth=1, color=self.temp_plot_colour_dict[temp], markersize=2,
            label=None
        )

    def _plot_adult_survival(self, i, j, sp, data_type):
        if sp in ['ad', 'ah']:
            # Plot the aduult_survival as just that, survival. So start with 100% and then decrease
            ax = self.axes[i][j]
            # We will need to do an ax.errorbar for each of the temperatures
            param_to_plot = "survival_percent"
            working_df = self.survival_df[
                (self.survival_df['adult_recruit'] == 'adult') & (self.survival_df['species'] == sp)]
            working_df['survival_percent'] = (working_df['survival'] / 5) * 100
            working_df = working_df[working_df[param_to_plot].notnull()]
            if working_df.index.to_list():
                self._plot_a_set_of_line_data(ax, data_type, sp, working_df, param_to_plot)
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
        xl_df = pd.read_excel(io=self.physiological_data_path, sheet_name=f'{sp}_recruit_size_fv_fm_bh')
        for ind in xl_df.index:
            self._populate_fv_fm_sample_rows(ind=ind, fv_fm_size_dd=fv_fm_size_dd, xl_df=xl_df, param='yield', sp=sp)
            self._populate_fv_fm_sample_rows(ind, fv_fm_size_dd, xl_df, 'cylinder_vol', sp=sp)
        return fv_fm_size_dd

    def _populate_recruit_survival(self, sp, survival_data):
        xl_df = pd.read_excel(io=self.physiological_data_path, sheet_name=f'{sp}_recruit_survival_bh')
        self._pop_survival_data(survival_data=survival_data, xl_df=xl_df, species=sp,
                                adult_recruit='recruit', time_unit='months', exp_type='main')

    def _populate_adult_survival(self, sp, survival_data):
        xl_df = pd.read_excel(io=self.physiological_data_path, sheet_name=f'{sp}_adult_survival_bh')
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

    def _get_ident_from_df_no_time(self, ind, df):
        ident = f"{df.at[ind, 'temperature']}_{df.at[ind, 'tank']}_{df.at[ind, 'rack']}_" \
                f"{df.at[ind, 'rack_row']}_{df.at[ind, 'rack_col']}"
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

    def _create_colour_list(self,
                            sq_dist_cutoff=None, mix_col=None, num_cols=50, time_out_iterations=10000,
                            avoid_black_and_white=True):
        new_colours = []
        min_dist = []
        attempt = 0
        while len(new_colours) < num_cols:
            attempt += 1
            # Check to see if we have run out of iteration attempts to find a colour that fits into the colour space
            if attempt > time_out_iterations:
                sys.exit('Colour generation timed out. We have tried {} iterations of colour generation '
                         'and have not been able to find a colour that fits into your defined colour space.\n'
                         'Please lower the number of colours you are trying to find, '
                         'the minimum distance between them, or both.'.format(attempt))
            if mix_col:
                r = int((random.randint(0, 255) + mix_col[0]) / 2)
                g = int((random.randint(0, 255) + mix_col[1]) / 2)
                b = int((random.randint(0, 255) + mix_col[2]) / 2)
            else:
                r = random.randint(0, 255)
                g = random.randint(0, 255)
                b = random.randint(0, 255)

            # now check to see whether the new colour is within a given distance
            # if the avoids are true also
            good_dist = True
            if sq_dist_cutoff:
                dist_list = []
                for i in range(len(new_colours)):
                    distance = (new_colours[i][0] - r) ** 2 + (new_colours[i][1] - g) ** 2 + (
                            new_colours[i][2] - b) ** 2
                    dist_list.append(distance)
                    if distance < sq_dist_cutoff:
                        good_dist = False
                        break
                # now check against black and white
                d_to_black = (r - 0) ** 2 + (g - 0) ** 2 + (b - 0) ** 2
                d_to_white = (r - 255) ** 2 + (g - 255) ** 2 + (b - 255) ** 2
                if avoid_black_and_white:
                    if d_to_black < sq_dist_cutoff or d_to_white < sq_dist_cutoff:
                        good_dist = False
                if dist_list:
                    min_dist.append(min(dist_list))
            if good_dist:
                new_colours.append((r, g, b))
                attempt = 0

        return new_colours

    def _get_pre_def_colour_dict(self):
        """These are the top 40 most abundnant named sequences. I have hardcoded their color."""
        return {
            'A1': "#FFFF00", 'C3': "#1CE6FF", 'C15': "#FF34FF", 'A1bo': "#FF4A46", 'D1': "#008941",
            'C1': "#006FA6", 'C27': "#A30059", 'D4': "#FFDBE5", 'C3u': "#7A4900", 'C42.2': "#0000A6",
            'A1bp': "#63FFAC", 'C115': "#B79762", 'C1b': "#004D43", 'C1d': "#8FB0FF", 'A1c': "#997D87",
            'C66': "#5A0007", 'A1j': "#809693", 'B1': "#FEFFE6", 'A1k': "#1B4400", 'A4': "#4FC601",
            'A1h': "#3B5DFF", 'C50a': "#4A3B53", 'C39': "#FF2F80", 'C3dc': "#61615A", 'D4c': "#BA0900",
            'C3z': "#6B7900", 'C21': "#00C2A0", 'C116': "#FFAA92", 'A1cc': "#FF90C9", 'C72': "#B903AA",
            'C15cl': "#D16100", 'C31': "#DDEFFF", 'C15cw': "#000035", 'A1bv': "#7B4F4B", 'D6': "#A1C299",
            'A4m': "#300018", 'C42a': "#0AA6D8", 'C15cr': "#013349", 'C50l': "#00846F", 'C42g': "#372101"}

    def _get_colour_list(self):
        colour_list = [
            "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062",
            "#0CBD66",
            "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0",
            "#BEC459",
            "#456648", "#0086ED", "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9",
            "#FF913F",
            "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
            "#7900D7",
            "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F",
            "#201625",
            "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55", "#0089A3",
            "#CB7E98",
            "#A4E804", "#324E72", "#6A3A4C", "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489",
            "#806C66",
            "#222800", "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
            "#C895C5",
            "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58", "#7A7BFF", "#D68E01",
            "#353339",
            "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D", "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A",
            "#001325",
            "#02525F", "#0AA3F7", "#E98176", "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75",
            "#8D8546",
            "#9695C5", "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
            "#00005F",
            "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01", "#6B94AA", "#51A058",
            "#A45B02",
            "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966", "#64547B", "#97979E", "#006A66", "#391406",
            "#F4D749",
            "#0045D2", "#006C31", "#DDB6D0", "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE",
            "#C6DC99",
            "#203B3C", "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
            "#C6005A",
            "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183", "#AA5199", "#B5D6C3",
            "#A38469",
            "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433", "#789EC9", "#6D80BA", "#953F00", "#5EFF03",
            "#E4FFFC",
            "#1BE177", "#BCB1E5", "#76912F", "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213",
            "#A76F42",
            "#89412E", "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
            "#BDC9D2",
            "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00", "#061203", "#DFFB71",
            "#868E7E",
            "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66", "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C",
            "#00B57F",
            "#545C46", "#866097", "#365D25", "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"]
        return colour_list

schupp_figures().plot_fig_1()

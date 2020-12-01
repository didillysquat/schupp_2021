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
        self.data_types = ['adult_survival', 'adult_zooxs', 'recruit_survival', 'recruit_growth', 'recruit_fv_fm']

        # Dataframe
        # We want to create one single dataframe that holds all of the physiology data.
        # We will hold the zooxs data in seperate dataframes that we will create from the SP outputs
        # TODO refactor
        # For the physiology data we will need a dataframe for survival and one for growth and fv_fm
        # This is because survival was assesed as a batch (i.e. 3 dead out of 5), while the growth and fv_fm
        # were done on a persample basis.

        # The survival df will have species, adult_recruit, time_value, time_unit, tank, temperature, exp_type
        # NB exp type will be 'main' or 'delayed'.


        # The physiology dataframe will have species, adult_recruit, tank, rack, plate_row, plate_no
        self.phys_df = pd.DataFrame(columns=[
            'species', 'adult_recruit', 'time_value', 'time_unit', 'tank', 'rack',
            'temperature', 'exp_type', 'plate_col', 'plate_row', 'cyl_vol', 'fv_fm'
        ])

        # get the survial data
        survival_data = []
        for sp in self.species_short:
            for data_type in self.data_types:
                if data_type == 'adult_survival':
                    xl_df = pd.read_excel(io=self.physiological_data_path, sheet_name=f'{sp}_adult_survival_bh')
                    self._pop_survival_data(survival_data=survival_data, xl_df=xl_df, species=sp,
                                    adult_recruit='adult',time_unit='day', exp_type='main')

                elif data_type == 'recruit_survival':
                    xl_df = pd.read_excel(io=self.physiological_data_path, sheet_name=f'{sp}_recruit_survival_bh')
                    self._pop_survival_data(survival_data=survival_data, xl_df=xl_df, species=sp,
                                            adult_recruit='recruit', time_unit='month', exp_type='main')
                foo = 'bar'
        self.survival_df = pd.DataFrame(
            columns=['species', 'adult_recruit', 'time_value', 'time_unit', 'tank', 'temperature', 'exp_type',
                     'survival'])
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

schupp_figures()

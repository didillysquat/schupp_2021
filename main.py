#!/usr/bin/env python3
"""
All code relating to the Schupp paper.
First task is to get all of the sequencing files into a single datasheet for analysis in SymPortal.
Currently the sequencing files have been analysed in two separate analyses as the sequencing data was available
at two time points.
Sam has sent me over some big excel sheets that contain all of the information we need:
esentially linking sample names to sequencing files.
The sample names were previously generated from the seq_files. We will continue to do it this way.

The Big excel files that Sam sent over was not really usable. I sent him a different datasheet for him to fill out
that we will now work with. It is called all_schupp_data_sam_filled. We will use this to make a
SymPortal datasheet so that we can get the samples back into the analysis.
"""

import pandas as pd

class SchuppAnalysis:
    def __init__(self):
        self.file_list_one_path = '/Users/humebc/Google_Drive/projects/schupp_larval/file_list_one.txt'
        self.file_list_two_path = '/Users/humebc/Google_Drive/projects/schupp_larval/file_list_two.txt'
        self.seq_filename_list = self._make_filename_list()

        self.path_to_info_spread = '/Users/humebc/Google_Drive/projects/schupp_larval/all_schupp_data_sam_filled.xlsx'
        self.info_df = self._make_info_df()

        # Make a dataframe that we can output as a csv that can then be pasted into the Symportal df
        self.sp_df_headers = "sample_name	fastq_fwd_file_name	fastq_rev_file_name	sample_type	host_phylum	host_class	host_order	host_family	host_genus	host_species	collection_latitude	collection_longitude	collection_date	collection_depth".split()
        self.sp_df = self._make_sp_df()

    def _make_sp_df(self):
        # First check to make sure that all of the sample_names (i.e. the index of the info_df) are found in the
        # file lists
        ind_to_seq_dict = {}
        for index_name in self.info_df.index:
            match = [candidate for candidate in self.seq_filename_list if index_name in candidate]

            if len(match) != 2:
                raise RuntimeError(f'{index_name} does not match exactly 2 seq files')
            else:
                if 'R1' in match[0]:
                    ind_to_seq_dict[index_name] = match
                else:
                    ind_to_seq_dict[index_name] = [match[1], match[0]]

        # If we get here then we know we have exact matches for each of the sample names
        # Now populate a 2D list that will be convereted into a df

        df_list = []
        for index_name in self.info_df.index:
            row = []
            row.append(index_name)
            row.extend(ind_to_seq_dict[index_name])

            # Sample type
            if self.info_df.at[index_name, 'age'] != 'adult':
                row.append('coral_aquarium')
            elif self.info_df.at[index_name, 'age'] == 'adult':
                row.append('coral_field')

            # Taxa info
            if self.info_df.at[index_name, 'Species'].rstrip() == 'A. digitifera':
                row.extend(['Cnidaria', 'Anthozoa', 'Scleractinia', 'Acroporidae', 'Acropora', 'digitifera'])
            elif self.info_df.at[index_name, 'Species'].rstrip() == 'A. hyacinthus':
                row.extend(['Cnidaria', 'Anthozoa', 'Scleractinia', 'Acroporidae', 'Acropora', 'surculosa'])
            elif self.info_df.at[index_name, 'Species'].rstrip() == 'L. purpurea':
                row.extend(['Cnidaria', 'Anthozoa', 'Scleractinia', 'incertae sedis', 'Leptastrea', 'purpurea'])
            elif self.info_df.at[index_name, 'Species'].rstrip() == 'P. damicornis':
                row.extend(['Cnidaria', 'Anthozoa', 'Scleractinia', 'Pocilloporidae', 'Pocillopora', 'damicornis'])
            else:
                raise RuntimeError(f'Unrecognised species of {self.info_df.at[index_name, "Species"]} for {index_name}')

            # Additional meta info
            # So that we can use this df as a single meta info file when making the figures etc we will include
            # info on temp tank and age.
            row.append(self.info_df.at[index_name, 'temp'])
            row.append(self.info_df.at[index_name, 'tank'])
            row.append(self.info_df.at[index_name, 'age'])
            df_list.append(row)
        headers_with_meta = self.sp_df_headers[:-4] + ['temp', 'tank', 'age']
        sp_df = pd.DataFrame(df_list, columns=headers_with_meta)
        sp_df.to_csv('sp_df.csv', index=False)
        return sp_df

    def _make_info_df(self):
        info_df = pd.read_excel(
            io=self.path_to_info_spread, header=0)
        assert(len(info_df['name'].index) == len(set(info_df['name'].index)))
        info_df.set_index('name', inplace=True, drop=True)
        return info_df

    def _make_filename_list(self):
        """Read in the filenames from the directory listings and put them into a single list"""
        filename_list = []
        with open(self.file_list_one_path, 'r') as f:
            filename_list.extend([line.rstrip().split()[8].replace('*', '') for line in f])

        with open(self.file_list_two_path, 'r') as f:
            filename_list.extend([line.rstrip().split()[8].replace('*', '') for line in f])

        return filename_list

sa = SchuppAnalysis()
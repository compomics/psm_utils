"""Interface to MS Amanda CSV result files."""

import logging
import os
import re
from functools import cmp_to_key
from typing import Dict, List, Optional, Tuple, Type, Union

import click
import numpy as np
import pandas as pd

from ms2rescore.peptide_record import PeptideRecord
from ms2rescore._exceptions import ModificationParsingError

logger = logging.getLogger(__name__)


@pd.api.extensions.register_dataframe_accessor("msamanda")
class MSAmandaAccessor:
    """ Pandas extension for MS Amanda csv files """

    default_columns = {
        'Scan Number', 
        'Title', 
        'Sequence', 
        'Modifications',
        'Protein Accessions', 
        'Amanda Score', 
        'Weighted Probability', 
        'Rank',
        'm/z', 
        'Charge', 
        'RT', 
        'Nr of matched peaks', 
        'Filename'
    }

    # Constructor
    def __init__(self,pandas_object) -> None:
        """ Pandas extenstion for MS Amanda csv files """
        self._obj = pandas_object
        self.invalid_amino_acids = r"[BJOUXZ]"

        self.modification_mapping = {}
        self.fixed_modifications = {}

    @classmethod
    def from_file(
        cls,
        path_to_MSAmanda_file: Union[str,os.PathLike],

    ) -> pd.DataFrame:
        """  """

        msamanda_df = pd.read_csv(path_to_MSAmanda_file, sep="\t",comment="#")

        msamanda_df = msamanda_df.msamanda.filter_rank1_psms() 
        msamanda_df = msamanda_df.msamanda.remove_psms_with_invalid_amino_acids()

        return msamanda_df


    
    def filter_rank1_psms(self) -> pd.DataFrame:
        """ Filters the MS Amanda file for rank 1 PSMs """
        
        self._obj = self._obj[self._obj['Rank']==1]

        # To ensure only one hit per spectrum
        self._obj = self._obj.sort_values("Amanda Score", ascending=False)
        duplicate_indices = self._obj[self._obj.duplicated(["Scan Number"], keep="first")].index
        self._obj = self._obj.drop(duplicate_indices).sort_index().reset_index()

        return self._obj

    def remove_psms_with_invalid_amino_acids(self) -> pd.DataFrame:
        """ Removes PSMs that have invalid amino acids """

        invalid_indices = self._obj[self._obj["Sequence"].str.contains(self.invalid_amino_acids, regex=True)].index

        self._obj = self._obj.drop(index=invalid_indices).reset_index(drop=True)

        return self._obj
    
    @staticmethod
    def extract_peprec_mod_notation_elements(mod_tuple,modification_mapping):
        """ Extracts PEPREC modification notation elements of the tuple and returns them as a list """
        
        # Map MSAmanda name to MS²Rescore name
        if mod_tuple[2] == "Phospho":
            phospho_name = mod_tuple[2]+mod_tuple[0]
            mod_name = modification_mapping[phospho_name]
        else:
            mod_name = modification_mapping[mod_tuple[2]]

        # If the peptide has (a) modification(s) on the N- or C-terminal
        if mod_tuple[1] == "-term":
            if mod_tuple[0] == "N":
                return ["0", mod_name]
            elif mod_tuple[0] == "C":
                return ["-1", mod_name]
        # If it is any other residue
        else: 
            return[mod_tuple[1], mod_name]

    @staticmethod
    def get_peprec_mod_notation_row(mod_tuples_list,modification_mapping):
        """ Returns a the PEPREC modification notation for one psm """
        
        peprec_mod_elements=[]

        for x in mod_tuples_list:
            peprec_mod_elements.extend(MSAmandaAccessor.extract_peprec_mod_notation_elements(x,modification_mapping))

        return "|".join(peprec_mod_elements)
    

    def get_peprec_modifications(self,modification_mapping) -> List:
        """ Returns a list of all modification in PEPREC format for the MS Amanda csv file """

        self._obj['Modifications'] = self._obj['Modifications'].fillna('-')

        modifications = self._obj['Modifications'].str.findall(r'([A-Z])(-term|\d+)\(([A-Za-z]+)\|([0-9.]+)\|(variable|fixed)\);?')
        

        peprec_mod_notations=[]

        for row in modifications:
            if row: # If the list of the row is not empty 
                peprec_mod_notations.append(MSAmandaAccessor.get_peprec_mod_notation_row(row,modification_mapping))
            else:
                peprec_mod_notations.append('-')

        return peprec_mod_notations

    def get_target_decoy_label(self) -> List:

        """ 
        Returns a list same length as the number of rows in self.obj_ 
        indicating whether the psm is a target (1) or decoy (-1) hit 

        """
        protein_accessions = self._obj['Protein Accessions']

        labels = []

        for row in protein_accessions:
            if "REV_" in row: # REV_ indicate that it is a decoy hit
                labels.append(-1)
            else:
                labels.append(1)

        return labels
    
    def get_protein_list(self) -> List:
        """
        Returns a list same length as the number of rows in self.obj_ 
        indicating whether the psm is a target (1) or decoy (-1) hit

        """

        protein_accessions = self._obj['Protein Accessions']

        protein_accessions_new = []

        for row in protein_accessions:
            protein_accessions_new.append(row.replace(";","|||"))

        return protein_accessions_new


    
    def to_peprec(self, modification_mapping=None) -> PeptideRecord:
        """ Get the PeptideRecord from the MS Amanda csv file."""

        if modification_mapping:
            self.modification_mapping = modification_mapping

        msamanda_peprec_format = pd.DataFrame(
            columns=[
                'spec_id',
                'peptide',
                'modifications',
                'charge',
                'protein_list',
                'psm_score',
                'observed_retention_time',
                'Label'
            ]
        )

        # Filling data frame 
        msamanda_peprec_format['spec_id'] = self._obj['Title']
        msamanda_peprec_format['peptide'] = self._obj['Sequence'].str.upper()
        msamanda_peprec_format['modifications'] = self.get_peprec_modifications(modification_mapping)
        msamanda_peprec_format['charge'] = self._obj['Charge']
        msamanda_peprec_format['protein_list'] = self.get_protein_list()
        msamanda_peprec_format['psm_score'] = self._obj['Amanda Score']
        msamanda_peprec_format['observed_retention_time'] = self._obj['RT']
        msamanda_peprec_format['Label'] = self.get_target_decoy_label()

        
        return PeptideRecord.from_dataframe(msamanda_peprec_format)

    
    
    def get_search_engine_features(self):
        """ Gets the features from the MS Amanda output file for rescoring by Percolator """

        features = self._obj[[
            'Scan Number',
            'Amanda Score',
            'Weighted Probability',
            'Charge',
            'm/z',
            'RT',
            'Nr of matched peaks',
        ]].rename(columns={
            'Scan Number': 'ScanNumber',
            'Amanda Score': 'Score',
            'Weighted Probability':'WeightedProbability',
            'm/z':'MassOverCharge',
            'RT':'RetentionTime',
            'Nr of matched peaks':'NumberOfMatchedPeaks',
        })

        return features


@click.command()
@click.argument("input-msamanda")
@click.argument("output-peprec")
def main(**kwargs):
    """ Convert msamanda csv output file to PEPREC """
    msamanda_df = pd.DataFrame.msamanda.from_file(kwargs["input_psm_report"])
    peprec = msamanda_df.msamanda.to_peprec()
    peprec.to_csv(kwargs["output_peprec"])


if __name__ == "__main__":
    main()



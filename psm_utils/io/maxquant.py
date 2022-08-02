"""Interface to MaxQuant msms.txt PSM files."""

import csv
import logging
import os
import re
from functools import cmp_to_key
from itertools import compress
from pathlib import Path
from typing import Dict, List, Tuple, Union

import numpy as np
import pandas as pd

from psm_utils.exceptions import PSMUtilsException
from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PeptideSpectrumMatch
from psm_utils.psm_list import PSMList

logger = logging.getLogger(__name__)

MSMS_DEFAULT_COLUMNS = {
    "Raw file",
    "Scan number",
    "Charge",
    "Length",
    "Sequence",
    "Modified sequence",
    "Proteins",
    "Missed cleavages",
    "Mass",
    "Mass error [Da]",
    "Mass error [ppm]",
    "Reverse",
    "Retention time",
    "PEP",
    "Score",
    "Delta score",
    "Localization prob",
    "Matches",
    "Intensities",
    "Mass Deviations [Da]",
    "Mass Deviations [ppm]",
    "Intensity coverage",
    "id",
}


class MaxquantReader(ReaderBase):
    """Reader for MaxQuant msms.txt PSM files."""

    def __init__(
        self,
        filename: Union[str, Path],
    ) -> None:
        super().__init__(filename)

        self._rename_mapping = None
        self._mass_error_unit = None

        self._validate_msms()

    def __iter__(self):
        """Iterate over file and return PSMs one-by-one"""

        with open(self.filename) as msms_in:
            reader = csv.DictReader(msms_in, delimiter="\t")

            for psm_dict in reader:
                psm = self._get_peptide_spectrum_match(psm_dict)
                yield psm

    # TODO next does not work yet
    def __next__(self) -> PeptideSpectrumMatch:
        return super().__next__()

    def read_file(self) -> PSMList:
        """Read full MaxQuant msms.txt PSM file into a PSMList object."""
        msms_psms = []
        with open(self.filename) as msms_in:
            reader = csv.DictReader(msms_in, delimiter="\t")

            for psm_dict in reader:
                msms_psms.append(self._get_peptide_spectrum_match(psm_dict))

        return PSMList(msms_psms)

    def _validate_msms(self) -> None:
        with open(self.filename, "r") as msms_file:
            msms_reader = csv.DictReader(msms_file, delimiter="\t")
            self._evaluate_columns(msms_reader.fieldnames)
            self._set_mass_error_unit(msms_reader.fieldnames)
            self._rename_mapping = self._fix_column_case(msms_reader.fieldnames)

    @staticmethod
    def _evaluate_columns(columns) -> bool:
        """Case insensitive column evaluation msms file."""

        columns = list(map(lambda col: col.lower(), columns))
        column_check = [
            True if col.lower() in columns else False for col in MSMS_DEFAULT_COLUMNS
        ]
        if all(column_check):
            return True
        else:
            raise MsmsParsingError(
                f"Missing columns: {list(compress(columns, column_check))}"
            )

    @staticmethod
    def _fix_column_case(columns: List[str]) -> Dict[str, str]:
        """
        Create mapping for column names with the correct case.

        Using `_evaluate_columns`, we can load required columns in a case-insensitive
        manner. As a result, the column name case must be fixed for downstream usage.
        """
        case_mapping = {col.lower(): col for col in MSMS_DEFAULT_COLUMNS}
        required_col = list(map(lambda col: col.lower(), MSMS_DEFAULT_COLUMNS))
        rename_mapping = {
            case_mapping[col.lower()]: col
            for col in columns
            if col.lower() in required_col
        }
        return rename_mapping

    def _set_mass_error_unit(self, columns) -> None:
        """Get mass error unit from DataFrame columns."""
        if "Mass error [Da]" in columns:
            self._mass_error_unit = "Da"
        elif "Mass error [ppm]" in columns:
            self._mass_error_unit = "ppm"
        else:
            raise NotImplementedError(f"MSMS.txt mass error unit not supported.")

    def _get_peptide_spectrum_match(
        self, psm_dict: Dict[str, Union[str, float]]
    ) -> PeptideSpectrumMatch:
        """Return a PeptideSpectrumMatch object from maxquat msms PSM"""

        peptide = self._get_peptidoform(psm_dict["Modified sequence"])
        spectrum_id = psm_dict["Scan number"]
        run = psm_dict["Raw file"]
        is_decoy = psm_dict["Reverse"] == "+"
        score = float(psm_dict["Score"])
        precursor_charge = int(psm_dict["Charge"])
        precursor_mz = float(psm_dict["m/z"])
        retention_time = float(psm_dict["Retention time"])
        protein_list = psm_dict["Proteins"]
        source = "maxquant"
        provenance_data = ({"maxquant_filename": self.filename},)
        metadata = {
            "Delta Score": psm_dict["Delta score"],
            # "Missed cleavages": psm_dict["enzInt"],
            "Localization prob": psm_dict["Localization prob"],
            "PepLen": psm_dict["Length"],
            "Precursor Intensity": psm_dict["Precursor Intensity"],
            "dM": psm_dict[f"Mass error [{self._mass_error_unit}]"],
            "Matches": psm_dict["Matches"],
            "Intensities": psm_dict["Intensities"],
            "Intensity coverage": psm_dict["Intensity coverage"],
            f"Mass Deviations [{self._mass_error_unit}]": psm_dict[
                self._rename_mapping[f"Mass Deviations [{self._mass_error_unit}]"]
            ],
        }
        return PeptideSpectrumMatch(
            peptide=peptide,
            spectrum_id=spectrum_id,
            run=run,
            is_decoy=is_decoy,
            score=score,
            precursor_charge=precursor_charge,
            precursor_mz=precursor_mz,
            retention_time=retention_time,
            protein_list=protein_list,
            source=source,
            provenance_data=provenance_data,
            metadata=metadata,
        )

    def _get_peptidoform(
        self,
        modified_peptideseq: str,
    ) -> Peptidoform:
        """Return a peptido form"""

        return Peptidoform(self._parse_maxquant_modification(modified_peptideseq))

    def _parse_maxquant_modification(self, modified_seq):
        """Parse modified maxquant seq to proforma sequence"""

        # pattern to match open and closed round brackets
        pattern = re.compile(r"\(((?:[^)(]+|\((?:[^)(]+|\([^)(]*\))*\))*)\)")
        proforma_seq = modified_seq.strip("_")
        seq_length = len(proforma_seq) + 1

        for match in pattern.finditer(modified_seq):
            # take match object string, remove leading and trailing bracket
            # and escape remaining brackets
            se_mod_string = match[0][1:-1].replace("(", "\(").replace(")", "\)")

            # if N-term mod
            if match.start() == 1:
                proforma_seq = re.sub(
                    f"\({se_mod_string}\)",
                    f"[{match[0][1:-1]}]-",
                    proforma_seq,
                )

            # if C-term mod
            elif match.end() == seq_length:
                proforma_seq = re.sub(
                    f"\({se_mod_string}\)",
                    f"-[{match[0][1:-1]}]",
                    proforma_seq,
                )

            # if modification on amino acid
            else:
                proforma_seq = re.sub(
                    f"\({se_mod_string}\)",
                    f"[{match[0][1:-1]}]",
                    proforma_seq,
                )

        return proforma_seq


class MaxquantWriter(WriterBase):
    """Writer for MaxQuant msms.txt files (not implemented)."""

    def __init__(self, filename):
        raise NotImplementedError()

    def write_psm(self, psm: PeptideSpectrumMatch):
        """Write a single PSM to the MaxQuant msms.txt PSM file."""
        raise NotImplementedError()

    def write_file(self, psm_list: PSMList):
        """Write entire PSMList to the MaxQuant msms.txt PSM file."""
        raise NotImplementedError()


class ModificationParsingException(PSMUtilsException):
    """Modification parsing error."""


class MsmsParsingError(PSMUtilsException):
    """Msms parsing error"""


@pd.api.extensions.register_dataframe_accessor("msms")
class MSMSAccessor:
    """Pandas extension for MaxQuant msms.txt files."""

    _mass_error_unit = None

    def __init__(self, pandas_obj) -> None:
        """Pandas extension for MaxQuant msms.txt files."""
        self._obj = pandas_obj
        self._set_mass_error_unit()
        self.invalid_amino_acids = r"[BJOUXZ]"

    @classmethod
    def _evaluate_columns(cls, column: str) -> bool:
        """Case insensitive column evaluation for Pandas.read_csv usecols argument."""
        return column.lower() in [col.lower() for col in cls.default_columns]

    @classmethod
    def _fix_column_case(cls, columns: List[str]) -> Dict[str, str]:
        """
        Create mapping for column names with the correct case.

        Using `_evaluate_columns`, we can load required columns in a case-insensitive
        manner. As a result, the column name case must be fixed for downstream usage.
        """
        case_mapping = {col.lower(): col for col in cls.default_columns}
        rename_mapping = {col: case_mapping[col.lower()] for col in columns}
        return rename_mapping

    @classmethod
    def from_file(
        cls,
        path_to_msms: Union[str, os.PathLike],
        filter_rank1_psms: bool = True,
        validate_amino_acids: bool = True,
    ) -> pd.DataFrame:
        """
        Read msms.txt from file.

        Parameters
        ----------
        path_to_msms : str, os.Pathlike
            path to msms.txt file
        filter_rank1_psms : bool, optional
            filter for rank 1 PSMs
        validate_amino_acids : bool, optional
            remove PSMs where the sequence includes an invalid amino acid; required for
            MS2PIP compatibility

        Returns
        -------
        msms : ms2rescore.maxquant.MSMS
            MSMS object (pandas.DataFrame with additional methods)
        """

        msms_df = pd.read_csv(path_to_msms, sep="\t", usecols=cls._evaluate_columns)
        msms_df.rename(columns=cls._fix_column_case(msms_df.columns), inplace=True)
        if filter_rank1_psms:
            msms_df = msms_df.msms.filter_rank1_psms()
        if validate_amino_acids:
            msms_df = msms_df.msms.remove_invalid_amino_acids()
        return msms_df

    def _set_mass_error_unit(self) -> None:
        """Get mass error unit from DataFrame columns."""
        if "Mass error [Da]" in self._obj.columns:
            self._mass_error_unit = "Da"
        elif "Mass error [ppm]" in self._obj.columns:
            self._mass_error_unit = "ppm"
        else:
            raise NotImplementedError(f"MSMS.txt mass error unit not supported.")

    def filter_rank1_psms(self) -> pd.DataFrame:
        """Filter MSMS for rank 1 PSMs."""
        self._obj = self._obj.sort_values("Score", ascending=False)
        duplicate_indices = self._obj[
            self._obj.duplicated(["Raw file", "Scan number"], keep="first")
        ].index
        self._obj = self._obj.drop(duplicate_indices).sort_index().reset_index()

        logger.debug(
            f"Found {len(self._obj)} rank 1 PSMs of which "
            f"{len(self._obj[self._obj['Reverse'] == '+']) / len(self._obj):.0%} are "
            "decoy hits."
        )
        if len(duplicate_indices) > 0:
            logger.warning("Removed %i non-rank 1 PSMs.", len(duplicate_indices))

        return self._obj

    def remove_invalid_amino_acids(self) -> pd.DataFrame:
        """Remove invalid amino acids from MSMS."""
        invalid_indices = self._obj[
            self._obj["Sequence"].str.contains(self.invalid_amino_acids, regex=True)
        ].index
        self._obj = self._obj.drop(index=invalid_indices).reset_index(drop=True)

        if len(invalid_indices) > 0:
            logger.warning(
                "Removed %i PSMs with invalid amino acids.", len(invalid_indices)
            )

        return self._obj

    def _get_spec_id(self) -> pd.Series:
        """Get PEPREC-style spec_id."""
        return (
            self._obj["Raw file"]
            + "."
            + self._obj["Scan number"].astype(str)
            + "."
            + self._obj["Scan number"].astype(str)
        ).rename("spec_id")

    @staticmethod
    def _minus_one_compare_fn(mod_1, mod_2):
        """Custom comparision function where `-1` is always larger."""
        location_1 = mod_1[0]
        location_2 = mod_2[0]

        if location_1 == -1:
            if location_2 == -1:
                return 0
            else:
                return 1
        elif location_2 == -1:
            return -1
        else:
            return location_1 - location_2

    @staticmethod
    def _find_mods_recursively(mod_seq, pattern_mapping, regex_pattern, mod_list=None):
        """
        Find modifications in MaxQuant modified sequence recursively.

        Parameters
        ----------
        mod_seq : string
            MaxQuant modified sequence stripped of flanking amino acids and
            underscores
        pattern_mapping : dict[str, tuple]
            Mapping of modification pattern to name (e.g. `"(ox)": "Oxidation"`)
        regex_pattern : re.Pattern
            Compiled regex pattern containing all modification labels, including
            amino acid prefix (if not N/C-terminal)
        mod_list : list, optional
            List with modification positions and labels to recursively extend.

        """
        if not mod_list:
            mod_list = []

        # Recursively find matches
        match = re.search(regex_pattern, mod_seq)

        if match:
            pattern = match.group(0)
            mod_name = pattern_mapping[pattern]

            # Handle N/C-terminal modification locations
            if match.start() == 0:
                mod_location = 0
            elif match.end() == len(mod_seq):
                mod_location = -1
            else:
                mod_location = match.start()

            mod_list.append((mod_location, mod_name))

            # Remove current modification and recurse
            mod_seq = re.sub(regex_pattern, "", mod_seq, count=1)
            mod_list = MSMSAccessor._find_mods_recursively(
                mod_seq, pattern_mapping, regex_pattern, mod_list
            )

        # Validate that all modifications are found
        else:
            if not re.fullmatch(r"[A-Z]+", mod_seq):
                raise ModificationParsingException(
                    f"Coud not match remaining modification labels in sequence "
                    f"`{mod_seq}`. Ensure that all modifications are "
                    "configured in the MaxQuant `modification_mapping` setting."
                )

        return mod_list

    @staticmethod
    def _get_single_peprec_modification(
        sequence, modified_sequence, modification_mapping, fixed_modifications
    ):
        """
        Get single PEPREC-style modifications from MaxQuant modified sequence.
        """

        # Prepare modifications regex pattern
        pattern_mapping = {}
        for label, name in modification_mapping.items():
            pattern_mapping[f"({label})"] = name
        regex_pattern = re.compile(
            "|".join([re.escape(p) for p in pattern_mapping.keys()])
        )

        # Find variable modifications
        mod_list = MSMSAccessor._find_mods_recursively(
            modified_sequence, pattern_mapping, regex_pattern
        )

        # Add fixed modifications
        for aa, name in fixed_modifications.items():
            mod_list.extend([(m.start() + 1, name) for m in re.finditer(aa, sequence)])

        # Sort and format mod_list
        if mod_list:
            mod_string = "|".join(
                [
                    "|".join([str(x) for x in mod])
                    for mod in sorted(
                        mod_list, key=cmp_to_key(MSMSAccessor._minus_one_compare_fn)
                    )
                ]
            )
        else:
            mod_string = "-"

        return mod_string

    def get_peprec_modifications(
        self, modification_mapping=None, fixed_modifications=None
    ) -> List:
        """
        Get PEPREC-formatted modifications for full MSMS.

        Parameters
        ----------
        modification_mapping: dict
            Mapping used to convert the MaxQuant modification labels to
            PSI-MS modification names (e.g.
            `{"M", "Oxidation (M)"): "Oxidation"}`)
        fixed_modifications: dict
            Dictionary (`{aa: mod}`) with fixed modifications to be added to the
            peprec. E.g. `{'C': 'Carbamidomethyl'}`. MaxQuant output does not
            include modifications that were set as fixed during the search. The
            first tuple element contains the one-letter amino acid code. The
            second tuple element contains the full modification name, as listed
            in the MS²PIP configuration.
        """

        if not modification_mapping:
            modification_mapping = {}
        if not fixed_modifications:
            fixed_modifications = {}

        # Remove surrounding underscores
        if "_" in self._obj["Modified sequence"].iloc[0]:
            mod_sequences = self._obj["Modified sequence"].str.extract(
                "_(.*)_", expand=False
            )
        else:
            mod_sequences = self._obj["Modified sequence"]

        # Apply over PSMs
        peprec_mods = []
        for seq, mod_seq in zip(
            self._obj["Sequence"].to_list(), mod_sequences.to_list()
        ):
            peprec_mods.append(
                self._get_single_peprec_modification(
                    seq, mod_seq, modification_mapping, fixed_modifications
                )
            )
        return peprec_mods

    @staticmethod
    def _calculate_top7_peak_features(
        intensities: List, mass_errors: List
    ) -> Tuple[np.ndarray]:
        """
        Calculate "top 7 peak"-related search engine features.

        The following features are calculated:
        - mean_error_top7: Mean of mass errors of the seven fragment ion peaks with the
          highest intensities
        - sq_mean_error_top7: Squared MeanErrorTop7
        - stdev_error_top7: Standard deviation of mass errors of the seven fragment ion
          peaks with the highest intensities

        """
        if not (isinstance(intensities, list) and isinstance(mass_errors, list)):
            return np.nan, np.nan, np.nan

        else:
            intensities = [float(i) for i in intensities]
            mass_errors = [float(i) for i in mass_errors]

            indices_most_intens = np.array(intensities).argsort()[-1:-8:-1]
            mass_errors_top7 = [(mass_errors[i]) for i in indices_most_intens]
            mean_error_top7 = np.mean(mass_errors_top7)
            sq_mean_error_top7 = mean_error_top7**2
            stdev_error_top7 = np.std(mass_errors_top7)

            return mean_error_top7, sq_mean_error_top7, stdev_error_top7

    @staticmethod
    def _calculate_ion_current_features(
        matches: List, intensities: List, intensity_coverage: List
    ) -> Tuple[np.ndarray]:
        """
        Calculate ion current related search engine features.

        The following features are calculated:
        - ln_explained_ion_current: Summed intensity of identified fragment ions,
          divided by that of all fragment ions, logged
        - ln_nterm_ion_current_ratio: Summed intensity of identified N-terminal
          fragments, divided by that of all identified fragments, logged
        - ln_cterm_ion_current_ratio: Summed intensity of identified N-terminal
          fragments, divided by that of all identified fragments, logged
        - ln_ms2_ion_current: Summed intensity of all observed fragment ions, logged

        """
        pseudo_count = 0.00001
        if not isinstance(intensities, list):
            return np.nan, np.nan, np.nan, np.nan
        else:
            ln_explained_ion_current = intensity_coverage + pseudo_count
            summed_intensities = sum([float(i) for i in intensities])

            # Calculate ratio between matched b- and y-ion intensities
            y_ion_int = sum(
                [
                    float(intensities[i])
                    for i, m in enumerate(matches)
                    if m.startswith("y")
                ]
            )
            y_int_ratio = y_ion_int / summed_intensities

            ln_nterm_ion_current_ratio = (
                y_int_ratio + pseudo_count
            ) * ln_explained_ion_current
            ln_cterm_ion_current_ratio = (
                1 - y_int_ratio + pseudo_count
            ) * ln_explained_ion_current
            ln_ms2_ion_current = summed_intensities / ln_explained_ion_current

            out = [
                ln_explained_ion_current,
                ln_nterm_ion_current_ratio,
                ln_cterm_ion_current_ratio,
                ln_ms2_ion_current,
            ]

        return tuple([np.log(x) for x in out])

    def to_peprec(
        self,
        modification_mapping=None,
        fixed_modifications=None,
    ) -> pd.DataFrame:
        """
        Get PeptideRecord from MaxQuant msms.txt file.

        Parameters
        ----------
        modification_mapping: dict
            Mapping used to convert the two-letter MaxQuant modification labels to
            PSI-MS modification names.
        fixed_modifications: dict
            Dictionary ({aa: mod}) can contain fixed modifications to be added to the
            peprec. E.g. `{'C': 'Carbamidomethyl'}`, as the MaxQuant output does not
            include modifications that were set as fixed during the search. The first
            tuple element contains the one-letter amino acid code. The second tuple
            element contains the full modification name, as listed in the values of
            `modification_mapping`.

        """
        peprec = pd.DataFrame(
            columns=[
                "spec_id",
                "peptide",
                "modifications",
                "charge",
                "protein_list",
                "psm_score",
                "observed_retention_time",
                "Label",
                "Raw file",
            ]
        )
        peprec["spec_id"] = self._get_spec_id()
        peprec["peptide"] = self._obj["Sequence"]
        peprec["modifications"] = self.get_peprec_modifications(
            modification_mapping, fixed_modifications
        )
        peprec["charge"] = self._obj["Charge"]

        # Fill NaN values in Proteins column for decoy PSMs
        # But first check that NaN Proteins only occur for decoy PSMs, if so:
        # fill these without the "REV_"
        peprec["protein_list"] = self._obj["Proteins"].str.split(";")
        if (peprec["protein_list"].isna() & self._obj["Reverse"].isna()).any():
            req_cols = zip(
                self._obj["Proteins"],
                self._obj["Reverse"],
                self._obj["Modified sequence"],
            )
            peprec["protein_list"] = [
                [modseq] if (type(rev) == float) & (type(prot) == float) else prot
                for prot, rev, modseq in req_cols
            ]
        peprec["protein_list"] = peprec["protein_list"].fillna(
            "REV_" + self._obj["Modified sequence"]
        )

        peprec["psm_score"] = self._obj["Score"]
        peprec["observed_retention_time"] = self._obj["Retention time"]
        peprec["Label"] = self._obj["Reverse"].isna().apply(lambda x: 1 if x else -1)
        peprec["Raw file"] = self._obj["Raw file"]
        peprec.sort_values("spec_id", inplace=True)
        peprec.reset_index(drop=True, inplace=True)

        return peprec

    def get_search_engine_features(self):
        """
        Get search engine features from MSMS for Percolator rescoring.

        Percolator features are derived from the MSGF2PIN script. See table 1 of
        Percolator-MSGF+ article (doi.org/10.1021/pr400937n).
        """
        logger.debug("Calculating search engine features...")

        spec_id = self._get_spec_id()
        charge = self._obj["Charge"].rename("charge")

        directly_copied = self._obj[
            [
                "Score",
                "Delta score",
                "Localization prob",
                "Charge",
                "Mass",
                "Length",
                f"Mass error [{self._mass_error_unit}]",
                "Missed cleavages",
            ]
        ].rename(
            columns={
                "Score": "RawScore",
                "Delta score": "RawDeltaScore",
                "Localization prob": "RawModLocProb",
                "Length": "PepLen",
                f"Mass error [{self._mass_error_unit}]": "dM",
                "Charge": "ChargeN",
                "Missed cleavages": "enzInt",
            }
        )

        absdM = self._obj[f"Mass error [{self._mass_error_unit}]"].abs().rename("absdM")

        charges_encoded = pd.get_dummies(
            self._obj["Charge"], prefix="Charge", prefix_sep=""
        )

        top7_features = pd.DataFrame(
            [
                self._calculate_top7_peak_features(i, md)
                for i, md in zip(
                    self._obj["Intensities"].str.split(";"),
                    self._obj["Mass Deviations [Da]"].str.split(";"),
                )
            ],
            columns=["MeanErrorTop7", "sqMeanErrorTop7", "StdevErrorTop7"],
        )

        ion_current_features = pd.DataFrame(
            [
                self._calculate_ion_current_features(m, i, ic)
                for m, i, ic in zip(
                    self._obj["Matches"].str.split(";"),
                    self._obj["Intensities"].str.split(";"),
                    self._obj["Intensity coverage"],
                )
            ],
            columns=[
                "lnExplainedIonCurrent",
                "lnNTermIonCurrentRatio",
                "lnCTermIonCurrentRatio",
                "lnMS2IonCurrent",
            ],
        )

        features = (
            pd.concat(
                [
                    spec_id,
                    charge,
                    directly_copied,
                    absdM,
                    charges_encoded,
                    top7_features,
                    ion_current_features,
                ],
                axis=1,
            )
            .sort_values("spec_id")
            .reset_index(drop=True)
        )

        return features

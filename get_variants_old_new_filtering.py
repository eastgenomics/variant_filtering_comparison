import argparse
import concurrent.futures
import dxpy as dx
import numpy as np
import pandas as pd
import subprocess

from collections import defaultdict
from mergedeep import merge
from pathlib import Path


def parse_args():
    """
    Parse the command line arguments inputs given

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(
        description='Information necessary for optimised filtering'
    )

    parser.add_argument(
        '--cases',
        type=str,
        required=True,
        help="Path to the TSV of all of the reported cases"
    )

    parser.add_argument(
        '--dx_old',
        type=str,
        required=True,
        help=(
            'path to relevant folder(s) with old filtering in the DNAnexus '
            'test project'
        )
    )

    parser.add_argument(
        '--dx_new',
        type=str,
        required=True,
        help=(
            'path to relevant folder(s) with new filtering in the DNAnexus '
            'test project'
        )
    )

    parser.add_argument(
        '--old_path',
        type=str,
        required=True,
        help=(
            'Path on your local machine where the old filtering files will be '
            'downloaded to'
        )
    )

    parser.add_argument(
        '--new_path',
        type=str,
        required=True,
        help=(
            'Path on your local machine where the new filtering files will be '
            'downloaded to'
        )
    )

    parser.add_argument(
        '--out',
        type=str,
        required=True,
        help='Name of output spreadsheet of cases and filtering results'
    )

    args = parser.parse_args()

    return args


def find_files_in_project(search_folder):
    """
    Find all of the filtered VCF files in the relevant DNAnexus folder(s)
    in the testing project

    Parameters
    ----------
    search_folder : str
        the folder to search for VCFs in (if multiple, separated by commas)

    Returns
    -------
    all_vcfs : list
        list of dicts, each dict containing info about one VCF file
    """
    all_vcfs = []
    folder_paths = search_folder.split(',')

    # Search within each folder and add everything to list
    for folder in folder_paths:
        filtered_vcfs = list(dx.find_data_objects(
            project='project-GZ6g3BQ45B5j8YPb8QB8X5kF',
            folder=folder,
            name="*.vcf.gz",
            name_mode='glob',
            classname='file',
            describe=True
        ))

        all_vcfs.extend(filtered_vcfs)

        print(f"Found {len(filtered_vcfs)} VCFs in folder {folder}")

    print(f"Found {len(all_vcfs)} VCFs total in all folders")

    return all_vcfs


def make_vcf_dict(filtered_vcfs, vcf_filter_type):
    """
    Make a dictionary, where each sample is the key and the

    Parameters
    ----------
    filtered_vcfs : list
        list of dicts, each dict containing info about one VCF file
    vcf_filter_type : str
        'old' or 'new' depending on whether the VCFs have old or new filtering
        applied

    Returns
    -------
    vcf_file_dict : dict
        dict with sample as key and file info as nested dict
    Example:
    {
        'X123456': {
            'R58.4_Adult onset neurodegenerative disorder_P': {
                'old': [
                    {
                        'file_id': 'file-XYZ',
                        'name': 'X123456-GM123456-TWE-N-EGG4_markdup_recalibrated_Haplotyper_annotated.optimised_filtered.vcf.gz',
                        'name_with_panel': 'X123456-GM123456-TWE-N-EGG4_markdup_recalibrated_Haplotyper_annotated_R58.4.optimised_filtered.vcf.gz'
                    }
                ]
            }
        },
        'X234567'...
    }
    """
    vcf_file_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    # For each file, get the panel string it was tested for and add file
    # info to the dictionary under the relevant sample and panel
    for idx, file in enumerate(filtered_vcfs):
        print(f"Getting info for file {idx} of {len(filtered_vcfs)}")
        # Get file name and sample name from that
        name = file['describe']['name']
        exome_number = name.split('-')[0]

        # Get the panel the sample was tested for via the job input
        job_id = file['describe']['createdBy']['job']
        panel = dx.bindings.dxjob.DXJob(
            dxid=job_id).describe()['runInput']['panel_string']

        # Get the R code of the panel to create an updated filename for the VCF
        # with the panel included. This means later when we download we won't
        # just overwrite the same file for each sample if multiple exist
        r_code = panel.split('_')[0]
        file_base = name.removesuffix('.optimised_filtered.vcf.gz')
        updated_filename = f"{file_base}_{r_code}.optimised_filtered.vcf.gz"

        # Add to our dictionary
        vcf_file_dict[exome_number][panel][vcf_filter_type].append({
            'file_id': file['id'],
            'name': name,
            'name_with_panel': updated_filename
        })

    print(
        f"Found {len(vcf_file_dict.keys())} unique samples from those folders"
    )

    return vcf_file_dict


def create_file_dict(old_filter_path, new_filter_path):
    """
    Create a single dict containing info for each sample, listing VCF
    files for the sample and panel for both old and new filtering
    Parameters
    ----------
    old_filter_path : str
        name of folder(s) in DNAnexus for the old filter VCFs
    new_filter_path : str
        name of folder(s) in DNAnexus for the new filter VCFs

    Returns
    -------
    old_vcf_dict : dict
        dict with sample as key and file info for old + new filtering
        as nested dict
    Example (with file duplicates):
    {
        'X123456': {
            'R58.4_Adult onset neurodegenerative disorder_P': {
                'old': [
                    {
                        'file_id': 'file-XYZ',
                        'name': 'X123456-GM123456-TWE-N-EGG4_markdup_recalibrated_Haplotyper_annotated.optimised_filtered.vcf.gz',
                        'name_with_panel': 'X123456-GM123456-TWE-N-EGG4_markdup_recalibrated_Haplotyper_annotated_R58.4.optimised_filtered.vcf.gz'
                    },
                    {
                        'file_id': 'file-CDE',
                        'name': 'X123456-GM123456-TWE-N-EGG4_markdup_recalibrated_Haplotyper_annotated.optimised_filtered.vcf.gz',
                        'name_with_panel': 'X123456-GM123456-TWE-N-EGG4_markdup_recalibrated_Haplotyper_annotated_R58.4.optimised_filtered.vcf.gz'
                    }
                ],
                'new': [
                    {
                        'file_id': 'file-ABC',
                        'name': 'X123456-GM123456-TWE-N-EGG4_markdup_recalibrated_Haplotyper_annotated.optimised_filtered.vcf.gz',
                        'name_with_panel': 'X123456-GM123456-TWE-N-EGG4_markdup_recalibrated_Haplotyper_annotated_R58.4.optimised_filtered.vcf.gz'
                    },
                    {
                        'file_id': 'file-QRS',
                        'name': 'X123456-GM123456-TWE-N-EGG4_markdup_recalibrated_Haplotyper_annotated.optimised_filtered.vcf.gz',
                        'name_with_panel': 'X123456-GM123456-TWE-N-EGG4_markdup_recalibrated_Haplotyper_annotated_R58.4.optimised_filtered.vcf.gz'
                    }
                ]
            }
        },
        'X234567'...
    }
    """
    # Find all the files from the folder(s) containing the VCFs with old
    # filtering and then new filtering
    old_vcfs = find_files_in_project(old_filter_path)
    new_vcfs = find_files_in_project(new_filter_path)

    # Make a dict from each and then merge the dicts
    old_vcf_dict = make_vcf_dict(old_vcfs, 'old')
    new_vcf_dict = make_vcf_dict(new_vcfs, 'new')
    merge(old_vcf_dict, new_vcf_dict)

    return old_vcf_dict


def check_file_duplicates(vcf_dict):
    """
    Print any samples which have duplicates of files for the same clinical
    indication, only keep one set of old and new files per clinical indication
    per sample

    Parameters
    ----------
    vcf_dict : dict
        dict with sample as key, with nested dicts for each clinical indication
        per sample containing lists of files with old and new filtering

    Returns
    -------
    updated_dict : dict
        dict with sample as key and nested dicts for each clinical indication
        for that sample (a single dict per old and new filtering)
    Example:
    {
        'X123456': {
            'R58.4_Adult onset neurodegenerative disorder_P': {
                'old': {
                    'file_id': 'file-XYZ',
                    'name': 'X123456-GM123456-TWE-N-EGG4_markdup_recalibrated_Haplotyper_annotated.optimised_filtered.vcf.gz',
                    'name_with_panel': 'X123456-GM123456-TWE-N-EGG4_markdup_recalibrated_Haplotyper_annotated_R58.4.optimised_filtered.vcf.gz'
                },
                'new': {
                    'file_id': 'file-ABC',
                    'name': 'X123456-GM123456-TWE-N-EGG4_markdup_recalibrated_Haplotyper_annotated.optimised_filtered.vcf.gz',
                    'name_with_panel': 'X123456-GM123456-TWE-N-EGG4_markdup_recalibrated_Haplotyper_annotated_R58.4.optimised_filtered.vcf.gz'
                }
            }
        },
        'X234567'...
    }
    """
    updated_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))

    for sample, sample_info in vcf_dict.items():
        for clin_ind in sample_info:
            # Check how many 'old' filtering files there are for that clin ind
            no_of_old_files = len(sample_info[clin_ind]['old'])
            if no_of_old_files > 1:
                print(
                    f"Warning - {sample} has {no_of_old_files} VCFs for old "
                    f"filtering found for clinical indication {clin_ind}"
                )

            updated_dict[sample][clin_ind]['old'] = sample_info[clin_ind]['old'][0]

            # Check how many 'new' filtering files there are for that clin ind
            no_of_new_files = len(sample_info[clin_ind]['new'])
            if no_of_new_files > 1:
                print(
                    f"Warning - {sample} has {no_of_new_files} VCFs for new "
                    f"filtering found for clinical_indication {clin_ind}"
                )

            updated_dict[sample][clin_ind]['new'] = sample_info[clin_ind]['new'][0]

    return updated_dict


def download_files(old_path, new_path, sample_info):
    """
    If each file isn't already downloaded in the specified directory,
    download there

    Parameters
    ----------
    old_path : str
        full path to where the VCFs with old filtering will be downloaded to
    new_path : str
        full path to where the VCFs with new filtering will be downloaded to
    sample_info : dict
        the dictionary holding file info for a single sample
    """
    # Make dirs for old and new paths if they don't exist
    Path(old_path).mkdir(parents=True, exist_ok=True)
    Path(new_path).mkdir(parents=True, exist_ok=True)

    for clin_ind in sample_info:

        old_file_id = sample_info[clin_ind]['old']['file_id']
        old_file_name = sample_info[clin_ind]['old']['name_with_panel']
        file_path = Path(f"{old_path}/{old_file_name}")

        # If a file by that name has not already been downloaded in that dir
        # then download it using dx download, renaming it compared to DNAnexus
        # to include the panel name
        if not file_path.exists():
            dx.download_dxfile(old_file_id, f"{old_path}/{old_file_name}")

        # Do same for the new filtering files in a separate folder
        new_file_id = sample_info[clin_ind]['new']['file_id']
        new_file_name = sample_info[clin_ind]['old']['name_with_panel']
        file_path = Path(f"{new_path}/{new_file_name}")

        if not file_path.exists():
            dx.download_dxfile(new_file_id, f"{new_path}/{new_file_name}")


def concurrent_download(vcf_dict, workers, old_path, new_path):
    """
    Concurrently download VCF files to respective folders locally

    Parameters
    ----------
    vcf_dict : dict
        final dict with GM number as key and dict with single VCF and sample
        info as value
    workers : int
        number of workers
    old_path : str
        full path to where the VCFs with old filtering will be downloaded to
    new_path : str
        full path to where the VCFs with new filtering will be downloaded to
    """
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        concurrent_jobs = {
            executor.submit(
                download_files, old_path, new_path, sample_info
            ): sample_info for sample_info in vcf_dict.values()
        }
        for future in concurrent.futures.as_completed(concurrent_jobs):
            try:
                data = future.result()
            except Exception as exc:
                print(
                    "Error downloading data for "
                    f"{concurrent_jobs[future]}: {exc}"
                )


def get_PASS_variants(folder, name):
    """
    Get the HGVSc and HGVSp for the PASS variants (passed filtering)

    Parameters
    ----------
    folder : str
        name of the folder the file is in
    name : str
        name of the file

    Returns
    -------
    variants : str
        a string containing all of the PASS variants, each separated by a
        newline character so they display nicely later in an Excel spreadsheet
    """

    filepath = f"{folder}/{name}"

    variant_output = subprocess.run(
        f"bcftools query -i'FILTER=\"PASS\"' -f '%CSQ_HGVSc\t%CSQ_HGVSp\n' {filepath}",
        shell=True,
        capture_output=True
    )

    variants = variant_output.stdout.decode()

    return variants


def get_total_variants(folder, vcf_name):
    """
    Get the total number of variants in a VCF file

    Parameters
    ----------
    folder : str
        The full path of the folder containing the VCFs (old or new filtering)
    vcf_name : str
        The name of the VCF file

    Returns
    -------
    variants : int
        The number of variants in the VCF
    """
    filepath = f"{folder}/{vcf_name}"

    variant_output = subprocess.run(
        f"zgrep -v ^# {filepath} | wc -l",
        shell=True,
        capture_output=True
    )

    variants = variant_output.stdout.decode()

    return variants


def get_variant_info(vcf_dict, old_path, new_path):
    """
    Adds information about the total number of variants in the panel and the
    variants/variant counts with old and new filtering

    Parameters
    ----------
    vcf_dict : dict
        dict with sample as key and nested dicts for each clinical indication
        for that sample with a single dict per old and new filtering
    old_path : str
        path to folder holding all the VCFs with old filtering
    new_path : str
        path to folder holding all the VCFs with new filtering

    Returns
    -------
    vcf_dict : dict
        dict with sample as key and nested dicts for each clinical indication
        for that sample with added variants and counts
    Example:
    {
        'X123456': {
            'R58.4_Adult onset neurodegenerative disorder_P': {
                'total_variants': 514,
                'old': {
                    'file_id': 'file-XYZ',
                    'name': 'X123456-GM123456-TWE-N-EGG4_markdup_recalibrated_Haplotyper_annotated.optimised_filtered.vcf.gz',
                    'name_with_panel': 'X123456-GM123456-TWE-N-EGG4_markdup_recalibrated_Haplotyper_annotated_R58.4.optimised_filtered.vcf.gz',
                    'variants': 'NM_000165.5:c.*3dup\t.\nNM_004820.5:c.1478A>G\tNP_004811.1:p.Tyr493Cys\n',
                    'variant_count': 2
                },
                'new': {
                    'file_id': 'file-TUV',
                    'name': 'X123456-GM123456-TWE-N-EGG4_markdup_recalibrated_Haplotyper_annotated.optimised_filtered.vcf.gz',
                    'name_with_panel': 'X123456-GM123456-TWE-N-EGG4_markdup_recalibrated_Haplotyper_annotated_R58.4.optimised_filtered.vcf.gz'
                    'variants': 'NM_000165.5:c.*3dup\t.\n',
                    'variant_count': 1
                }
            }
        },
        'X234567'...
    }
    """
    for sample, sample_info in vcf_dict.items():
        for clin_ind in sample_info:
            # Get PASS variants with old filtering and add to dict
            old_filename = sample_info[clin_ind]['old']['name_with_panel']
            old_variants = get_PASS_variants(old_path, old_filename)
            sample_info[clin_ind]['old']['variants'] = old_variants
            # Add count from this to dict
            old_variant_list = old_variants.split('\n')
            old_var_list = list(filter(None, old_variant_list))
            sample_info[clin_ind]['old']['variant_count'] = len(old_var_list)

            # Get PASS variants with new filtering and add to dict
            new_filename = sample_info[clin_ind]['new']['name_with_panel']
            new_variants = get_PASS_variants(new_path, new_filename)
            sample_info[clin_ind]['new']['variants'] = new_variants
            # Add count from this to dict
            new_variant_list = new_variants.split('\n')
            new_var_list = list(filter(None, new_variant_list))
            sample_info[clin_ind]['new']['variant_count'] = len(new_var_list)

            total_variants = get_total_variants(new_path, new_filename)
            sample_info[clin_ind]['total_variants'] = int(total_variants.rstrip('\n'))

    return vcf_dict


def create_df(vcf_dict):
    """
    Create a df of the exome number, panel and variants so we can just merge
    this with the original df containing the reported variants

    Parameters
    ----------
    vcf_dict : dict
        dict of each sample and the routine and optimised variants and counts

    Returns
    -------
    variant_df : pd.DataFrame
        dataframe from the above dict
    """
    df_rows = []
    for sample, file_info in vcf_dict.items():
        for clin_ind in file_info:
            routine_variants = file_info[clin_ind]['old']['variants'].rstrip('\n')
            routine_count = file_info[clin_ind]['old']['variant_count']
            optimised_variants = file_info[clin_ind]['new']['variants'].rstrip('\n')
            optimised_count = file_info[clin_ind]['new']['variant_count']
            total_in_panel = file_info[clin_ind]['total_variants']

            df_rows.append({
                'ExomeNumber': sample,
                'CI': clin_ind,
                'Total_in_panel': total_in_panel,
                'Routine_variants': routine_variants,
                'Routine_count': routine_count,
                'Optimised_variants': optimised_variants,
                'Optimised_count': optimised_count,
                'Filtering_difference': routine_count - optimised_count
            })

    variant_df = pd.DataFrame(df_rows)

    return variant_df


def merge_cases_and_results(variant_df, cases_df):
    """
    Merge the variant and casees df to get one final df with reported
    variant(s) and routine and optimised filter variants

    Parameters
    ----------
    variant_df : pd.DataFrame
        dataframe with exome number, panel and the variants and counts
        for routine and optimised filtering
    cases_df : pd.DataFrame
        dataframe for each exome number and case with the originally
        reported variant and other info e.g. outcome for that case

    Returns
    -------
    merged_df : pd.DataFrame
        dataframe merged so for each case has the reported variant(s)
        for that case as well as the routine and optimised variants and
        counts
    """
    # Merge this with the variants we've found by exome no. and CI
    merged_df = variant_df.merge(
        cases_df, how='left', on=['ExomeNumber', 'CI']
    )

    merged_df.drop_duplicates(keep='first', inplace=True)

    return merged_df


def add_in_project_name_and_assay(variant_df):
    """
    Add in the name of the DX project the sample is from and the assay

    Parameters
    ----------
    variant_df : pd.DataFrame
        dataframe containing info about each case and variants found

    Returns
    -------
    reordered_df : pd.DataFrame
        dataframe with extra columns, re-ordered
    """

    variant_df['Project_name'] = variant_df.apply(
        lambda row: dx.describe(row['project_id'])['name'],
        axis=1
    )

    # Get assay from the project name and add new column
    variant_df['Assay'] = variant_df['Project_name'].str.split('_').str[-1]

    reordered_df = variant_df[[
        'ExomeNumber', 'Assay', 'HGVS1', 'Location1', 'Classification1',
        'HGVS2', 'Location2', 'Classification2', 'final_result', 'CI',
        'Short R-codes', 'Total_in_panel', 'Reported', 'Routine_variants',
        'Routine_count', 'Optimised_variants', 'Optimised_count',
        'Filtering_difference', 'project_id', 'file_id', 'Project_name'
    ]]

    return reordered_df


def add_number_of_reported_variants(variant_df):
    """
    Add count of reported variants based on content of HGVS1 and HGVS2 columns

    Parameters
    ----------
    variant_df : pd.DataFrame
        pandas dataframe with all the info about each case

    Returns
    -------
    variant_df : pd.DataFrame
        pandas dataframe with new column 'Reported' containing the number
        of variants reported
    """
    conditions = [
        ((variant_df['HGVS1'] == 'NMD') & (variant_df['HGVS2'] == 'NMD')),
        ((variant_df['HGVS1'] == 'NMD') & (variant_df['HGVS2'] != 'NMD')),
        ((variant_df['HGVS1'] != 'NMD') & (variant_df['HGVS2'] == 'NMD')),
        ((variant_df['HGVS1'] != 'NMD') & (variant_df['HGVS2'] != 'NMD'))
    ]

    values = [0, 1, 1, 2]

    variant_df['Reported'] = np.select(conditions, values)

    return variant_df


def main():
    args = parse_args()

    vcf_dict = create_file_dict(args.dx_old, args.dx_new)

    updated_vcf_dict = check_file_duplicates(vcf_dict)

    print("Downloading files")
    concurrent_download(
        updated_vcf_dict,
        8,
        args.old_path,
        args.new_path
    )

    print("Getting variants")
    vcf_dict_with_variants = get_variant_info(
        updated_vcf_dict,
        args.old_path,
        args.new_path
    )

    print("Creating dataframes")
    # Convert to df
    variant_df = create_df(vcf_dict_with_variants)
    cases_df = pd.read_csv(args.cases, sep='\t')
    # Merge the variant df with cases df which contains info on reported cases
    final_df = merge_cases_and_results(variant_df, cases_df)
    # Add the number of reported variants
    final_df = add_number_of_reported_variants(final_df)

    print("Getting project name and assay type")
    # Add name of the project the sample is from and the assay type
    final_df = add_in_project_name_and_assay(final_df)

    # Remove duplicates by ExomeNumber and CI
    final_df.drop_duplicates(
        subset=['ExomeNumber', 'CI'],
        keep='first',
        inplace=True
    )

    # Write out to CSV file
    final_df.to_excel(args.out, index=False)


if __name__ == "__main__":
    main()

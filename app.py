#!/usr/bin/env python3

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO import parse as seqio_parse
from collections import defaultdict
import gzip
from io import StringIO
import logging
import os
import pandas as pd
import plotly.express as px
from plotly import graph_objects as go
from scipy import cluster
import streamlit as st


###############################
# DEFINE ALL HELPER FUNCTIONS #
###############################
def read_data(
    dir_path
):
    """Read in data contained within a directory."""

    # Populate a dictionary with the output
    output = defaultdict(dict)

    # Walk the directory hierarchy
    for dirpath, dirnames, filenames in os.walk(dir_path):

        # Iterate over each file
        for filename in filenames:

            # Instead of checking file extensions, we'll
            # figure out the format based on the file contents

            # In contrast, the organism is determined by the extension
            organism, filetype, data = parse_file(
                os.path.join(dirpath, filename)
            )

            # If this filetype was not recognized
            if filetype is None:

                # Skip it
                continue
                
            # Otherwise
            else:

                logging.info(f"Read in {filetype} for {organism}")

                # Add it to the data
                output[organism][filetype] = data

    # Return the dict of all data
    return output


@st.cache
def parse_file(filepath, ignored_ext=('.zip', '.sh', '.DS_Store', ".docx")):
    """Parse the contents of a file."""

    # Make sure that the file exists
    assert os.path.exists(filepath), f"File not found: {filepath}"

    # Parse the organism name from the filepath
    org_name = parse_org_name(filepath)

    # If the file is any ignored type
    if filepath.endswith(ignored_ext):

        # Skip it and return none
        return org_name, None, None

    # If the file is gzip compressed
    if filepath.endswith(".gz"):
        logging.info(f"Reading {filepath} as gzip")

        # Open it with gzip
        with gzip.open(filepath, 'rt') as handle:

            # Try to
            try:

                # Read all of the contents
                lines = handle.readlines()

            # If the file is not unicode
            except UnicodeDecodeError:

                # Catch the error and return none
                return org_name, None, None

    # Otherwise
    else:
        logging.info(f"Reading {filepath}")

        # Open it as a text file
        with open(filepath, 'r') as handle:

            # Try to
            try:

                # Read all of the contents
                lines = handle.readlines()

            # If the file is not unicode
            except UnicodeDecodeError:

                # Catch the error and return none
                return org_name, None, None

    # Remove any empty header lines
    while len(lines) > 0:
        if len(lines[0]) <= 1:
            lines = lines[1:]
        else:
            break

    # If the first character is a '>'
    if lines[0][0] == ">":

        # Then it is a FASTA file
        return org_name, 'fasta', parse_fasta(lines)

    # If the first line starts with LOCUS
    elif lines[0].startswith("LOCUS"):

        # Then it is a GBK file
        return org_name, 'gbk', parse_gbk(lines)

    # If the first line contains a comma
    elif "," in lines[0]:

        # Parse the REBASE output
        return org_name, 'rebase', parse_rebase(lines)

    # Otherwise
    else:

        # The file type is not recognized
        return org_name, None, None    


def parse_fasta(lines):
    """Parse a FASTA file as a dict of {header: seq}."""

    return {
        header: seq
        for header, seq in SimpleFastaParser(
            StringIO("".join(lines))
        )
    }


def parse_gbk(lines):
    """Parse a GBK file as a DataFrame."""

    # Set up a list to fill with annotations
    data = []

    # Iterate over each record
    for record in seqio_parse(StringIO("".join(lines)), "genbank"):

        # Iterate over each feature
        for feature in record.features:

            data.append(
                dict(
                    record_id=record.id,
                    record_name=record.name,
                    record_description=record.description,
                    start=feature.location.start,
                    end=feature.location.end,
                    strand=feature.strand,
                    type=feature.type,
                    feature_id=feature.id,
                )
            )

            # Add all of the qualifiers for the feature
            for k, v in feature.qualifiers.items():
                data[-1][k] = v

    return pd.DataFrame(data)

def parse_rebase(lines):
    """Parse the epigenetic data from REBASE in CSV format."""

    return pd.read_csv(StringIO("".join(lines)))


def parse_org_name(filepath):
    """Parse the name of the organism from the name of the file."""

    # Strip away everything except the filename
    org_name = filepath.rsplit("/", 1)[-1]

    # If the file ends with ".gz"
    if filepath.endswith(".gz"):

        # Remove that extension
        org_name = org_name[:-3]

    # Iterate over a set of potential suffixes
    for suffix in [
        '.txt',
        '.gbk',
        '.gff',
        '.txt',
        '.fna',
        '.fasta',
        '.fa',
        '.csv'
    ]:
        # If the file ends with this suffix
        if org_name.endswith(suffix):

            # Remove it
            org_name = org_name[:-len(suffix)]

    # Return the organism name
    return org_name


@st.cache
def format_rebase_df(dirpath):

    # READ INPUT DATA FROM WORKING DIRECTORY
    data = read_data(dirpath)

    # Make a DataFrame with the organism name added
    df = pd.concat(
        [
            data[org_name]['rebase'].assign(
                organism=org_name
            )
            for org_name in data
            if 'rebase' in data[org_name]
        ]
    )

    # Format a label for the motif which combines multiple columns
    df = df.assign(
        motif_label=df.apply(format_motif_label, axis=1)
    )

    # Transform to a wide DataFrame
    df = df.pivot(
        index="organism",
        columns="motif_label",
        values="fraction"
    ).fillna(0)

    return df


def format_motif_label(r):
    """Format the label of each motif from the REBASE CSV file."""
    return f"{r.motifString}-pos{r.centerPos}-{'other' if r.modificationType == 'modified_base' else r.modificationType}"


def sort_table(df):
    """Sort the rows and columns of a table by linkage clustering."""

    return df.iloc[
        reordered_index(df),
        reordered_index(df.T)
    ]

def reordered_index(df, method="ward", metric="euclidean"):
    """Return the list of index positions following linkage clustering."""

    return cluster.hierarchy.leaves_list(
        cluster.hierarchy.linkage(
            df.values,
            method=method,
            metric=metric,
        )
    )


####################
# START THE SCRIPT #
####################

# FORMAT THE REBASE DATA
rebase_df = format_rebase_df(os.getcwd())

# TITLE OF THE APP
st.title('Epi-Pangenome Map')

# SIDEBAR MENUS
user_inputs = dict(
    hidden_genomes=st.sidebar.multiselect(
        "Hide Genomes",
        list(rebase_df.index.values)
    ),
    hidden_motifs=st.sidebar.multiselect(
        "Hide Motifs",
        list(rebase_df.columns.values)
    ),
    figure_height=st.sidebar.slider(
        "Figure Height",
        min_value=100,
        max_value=1200,
        step=1,
        value=600
    ),
    min_fraction=st.sidebar.slider(
        "Minimum Fraction per Motif",
        min_value=0.,
        max_value=1.,
        value=0.25,
        step=0.01,
        help="Only show motifs which are present in a single genome at this minimum threshold"
    ),
    min_prevalence=st.sidebar.slider(
        "Minimum Number of Genomes per Motif",
        min_value=1,
        max_value=rebase_df.shape[0],
        value=1,
        step=1,
        help="Only show motifs which are present in at least this many genomes"
    )
)

# MASK ANY SELECTED ROWS/COLUMNS
plot_df = rebase_df.drop(
    columns=user_inputs['hidden_motifs'],
    index=user_inputs['hidden_genomes'],
)

# MASK ANY MOTIFS WHICH DO NOT REACH THE MINIMUM THRESHOLD
plot_df = plot_df.drop(
    columns=[
        motif_name
        for motif_name, max_fraction in plot_df.max().items()
        if max_fraction < user_inputs['min_fraction']
    ]
)

# MASK ANY MOTIFS WHICH ARE NOT FOUND IN THE SPECIFIED NUMBER OF GENOMES
if user_inputs['min_prevalence'] > 1:

    plot_df = plot_df.reindex(
        columns=plot_df.columns.values[
            (plot_df > 0).sum() >= user_inputs['min_prevalence']
        ]
    )

# SORT THE ROWS/COLUMNS
plot_df = sort_table(plot_df)

# MAKE THE PLOT
st.write(
    go.Figure(
        go.Heatmap(
            x=plot_df.columns.values,
            y=plot_df.index.values,
            z=plot_df.values,
            colorscale="blues",
            text=plot_df.apply(
                lambda c: c.apply(
                    lambda v: f"Motif: {c.name}<br>Modification: {round(v * 100, 2):,}%"
                )
            ),
            hoverinfo="text"
        ),
        layout=dict(
            height=user_inputs['figure_height']
        )
    )
)
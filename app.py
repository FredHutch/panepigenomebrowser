#!/usr/bin/env python3

from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO import parse as seqio_parse
from collections import defaultdict
import gzip
from io import StringIO
from itertools import product
import logging
import os
import pandas as pd
import plotly.express as px
from plotly import graph_objects as go
from plotly.subplots import make_subplots
import re
from scipy import cluster
import streamlit as st
from zipfile import ZipFile


###############################
# DEFINE ALL HELPER FUNCTIONS #
###############################
@st.cache
def read_data(
    dir_path
):
    """Read in data contained within a directory."""

    # Populate a dictionary with the output
    output = defaultdict(dict)

    # By default, there are no additional annotations for the motifs or genomes
    addl_annots = []

    # Walk the directory hierarchy
    for dirpath, dirnames, filenames in os.walk(dir_path):

        # Iterate over each file
        for filename in filenames:

            logging.info(f"Considering {filename}")

            # If this is a CSV, then it may be a motif or annotation table
            if filename.endswith('.csv'):

                # Try to read it in
                try:
                    addl_annots.append(
                        pd.read_csv(
                            os.path.join(dirpath, filename),
                            index_col=0
                        )
                    )

                # If there are any errors parsing the file
                except:
                    # Just skip it
                    continue
            
            # If this is a zip file
            elif filename.endswith('.zip'):

                # Parse the contents
                for organism, filetype, data in parse_zip(
                    os.path.join(dirpath, filename)
                ):

                    # Add it to the data
                    output[organism][filetype] = data

            # Otherwise, if it is not a zip file
            else:

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

    # Get all of the annotations per-motif which agree across
    # all of the rebase files
    motif_annot = pd.concat([
        org_data['rebase'].drop(
            columns=[
                "org_name",
                "org_num",
                "percent_detected",
                "coverage",
                "comp_percent_detected",
                "pb_comment",
                "complete_genome",
                "enz_comment",
                "enz_name"
            ],
            errors="ignore"
        )
        for org, org_data in output.items()
        if 'rebase' in org_data
    ]).groupby(
        "rec_seq"
    ).apply(
        lambda d: pd.DataFrame([{
            col_name: col_values.dropna().unique()[0]
            for col_name, col_values in d.items()
            if col_values.dropna().unique().shape[0] == 1
        }])
    ).drop(
        columns=["rec_seq"]
    ).reset_index(
    ).drop(
        columns=["level_1"]
    ).assign(
        motif_length = lambda d: d.rec_seq.apply(len)
    ).set_index("rec_seq")

    # Now add in any additional annotations provided by the user
    for df in addl_annots:
        
        # Iterate over the columns of data provided by the user
        for col_name, col_values in df.items():

            # If any of the index labels from the annotation table
            # match the motifs used in the analysis
            if any([pd.notnull(col_values.get(i)) for i in motif_annot.index.values]):
            
                # Assign the values provided by the user to the rebase data
                motif_annot = motif_annot.assign(
                    **{
                        col_name: col_values
                    }
                )

    # Finally, reformat the index used to label the motif
    # so that it includes both the enzyme type and subtype
    # as indicated by the rebase data
    motif_annot = motif_annot.assign(
        type_label=motif_annot.apply(
            format_enzyme_type_label,
            axis=1
        )
    )

    # Make a DataFrame which can be used to populate annotations for the genomes
    genome_annot = pd.DataFrame(
        index=list(output.keys())
    )

    # Now add in any additional annotations provided by the user
    for df in addl_annots:
        
        # Iterate over the columns of data provided by the user
        for col_name, col_values in df.items():

            # If any of the index labels from the annotation table
            # match the motifs used in the analysis
            if any([pd.notnull(col_values.get(i)) for i in genome_annot.index.values]):
            
                # Assign the values provided by the user to the rebase data
                genome_annot = genome_annot.assign(
                    **{
                        col_name: col_values
                    }
                )

    # Return the dict of all data, along with the additional motif annotations, if any
    return output, motif_annot, genome_annot


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

    return parse_lines(org_name, lines)
    
@st.cache
def parse_zip(archivepath, ignored_ext=('.zip', '.sh', '.DS_Store', ".docx")):
    """Parse the contents of a zip archive."""

    logging.info(f"Parsing data from ZIP archive: {archivepath}")

    # Make sure that the file exists
    assert os.path.exists(archivepath), f"Archive not found: {archivepath}"

    # Set up a list of outputs
    outputs = []

    # Open the archive
    logging.info(f"Reading {archivepath} archive")
    with ZipFile(archivepath, 'r') as myzip:

        # Iterate over each file in the archive
        for filepath in myzip.namelist():
            logging.info(f"Parsing {filepath} from archive")

            # If this path is pointing to a folder inside the archive
            if filepath.endswith('/'):

                # Log it
                logging.info("This path appears to point to a folder, skipping")

                # Skip it
                continue

            # Parse the organism name from the filepath
            org_name = parse_org_name(filepath)

            # If the file is any ignored type
            if filepath.endswith(ignored_ext):

                # Skip it
                continue

            # If the file is gzip compressed
            if filepath.endswith(".gz"):

                # Open it with gzip
                with gzip.open(
                    myzip.open(filepath, 'r'),
                    'rt'
                ) as handle:

                    # Try to
                    try:

                        # Read all of the contents
                        lines = handle.readlines()

                    # If the file is not unicode
                    except UnicodeDecodeError:

                        # Catch the error and skip it
                        continue

            # Otherwise
            else:

                # Open it as a text file
                with myzip.open(filepath, 'r') as handle:

                    # Try to
                    try:

                        # Read all of the contents
                        lines = [
                            line.decode("utf-8").rstrip("\r")
                            for line in handle.readlines()
                        ]

                    # If the file is not unicode
                    except UnicodeDecodeError:

                        # Catch the error and skip it
                        continue

            # Parse the file
            org_name, filetype, dat = parse_lines(org_name, lines)

            # If the filetype was recognized
            if filetype is not None:

                # Add it to the list
                outputs.append((org_name, filetype, dat))

    # Return all of the data from the archive
    return outputs

def parse_lines(org_name, lines):
    """Parse the contents of a file, formatted as a list of lines."""

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

    # If the first line the <* delimiter
    elif lines[0].startswith("<*"):

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
                data[-1][k] = ", ".join(v) if isinstance(v, list) else v

    return pd.DataFrame(data)

def parse_rebase(lines):
    """Parse the epigenetic data from REBASE."""

    # Populate a list with dicts, each entry being an enzyme
    enzymes = list()
    enzyme = dict()

    # Keep track of the recognition sequences which have been detected thus far
    added_rec_seqs = set()

    # Iterate over each line
    for line in lines:

        # If the line is empty
        if len(line) <= 1:

            # Skip it
            continue

        # If the line has a top-level description field
        elif line.startswith("<*"):

            # Skip it
            continue

        # If the line is a field ending
        elif line.startswith("<>"):

            # If the current enzyme has content
            if len(enzyme) > 0:

                # If the recognition sequence is new
                if "rec_seq" in enzyme and enzyme["rec_seq"] not in added_rec_seqs:

                    # Add it to the list
                    enzymes.append(enzyme)

                    # Record that we've added this recognition sequence
                    added_rec_seqs.add(enzyme["rec_seq"])

            # Start a new blank entry for the next enzyme
            enzyme = dict()

        # If the line contains a key-value pair
        elif line.startswith("<"):

            # Make sure that there is a matching '>'
            assert '>' in line, f"Expected a '>' character in: {line}"

            # Parse the key and value
            key, value = line[1:].rstrip("\n").split(">", 1)

            # Add the key and value to the dict
            if value.isnumeric():
                enzyme[key] = int(value)
            elif value.isdecimal():
                enzyme[key] = float(value)
            else:
                enzyme[key] = value


    # At the end of reading all of the lines
    # If there is a field remaining
    if len(enzyme) > 0:

        # If the recognition sequence is new
        if "rec_seq" in enzyme and enzyme["rec_seq"] not in added_rec_seqs:

            # Add it to the list
            enzymes.append(enzyme)

    # Reformat the list of enzymes as a DataFrame
    return pd.DataFrame(enzymes)


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
    data, motif_annot, genome_annot = read_data(dirpath)

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
    
    # Only keep those records with the `percent_detected` field
    df = df.loc[
        df.percent_detected.notnull()
    ]

    # Add a single combined name to use for the display
    df = df.assign(
        text = df.apply(
            format_rebase_text,
            axis=1
        )
    )

    # Transform to a wide DataFrame containing floats for `percent_detected`
    percent_detected = df.pivot(
        index="organism",
        columns="rec_seq",
        values="percent_detected"
    ).fillna(
        0
    ).applymap(
        float
    )

    # Also make a wide table with the complete set of REBASE outputs
    text_df = df.pivot(
        index="organism",
        columns="rec_seq",
        values="text"
    ).fillna(
        ""
    ).reindex(
        index=percent_detected.index.values,
        columns=percent_detected.columns.values,
    )

    return percent_detected, text_df, motif_annot, genome_annot


def format_enzyme_type_label(r):
    """Format the label of each motif from the REBASE TXT file."""

    # If the enzyme type is present
    if pd.notnull(r.get("enz_type")):

        # And the subtype is present
        if pd.notnull(r.get("sub_type")):

            # Make a combined label
            return f"Type {int(r.enz_type)}{r.sub_type}"

        # Without the subtype
        else:
            return f"Type {int(r.enz_type)}"
    
    # Without any type information
    return None


def format_rebase_text(r):
    """Format the string which contains the complete set of REBASE information."""
    return "<br>".join(
        [
            f"{k}: {v}"
            for k, v in r.items()
            if k not in ['organism', 'enzyme_name']
        ]
    )


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


def format_display(plot_value_df, plot_text_df, motif_annot, genome_annot, user_inputs):
    """Format the display."""

    # If the user wants to sort the enzymes by an annotation
    if user_inputs["sort_motifs_by"] == "Motif Annotations":

        assert len(user_inputs["annot_motifs_by"]) > 0, "Must specify motif annotations for sorting"

        # Sort the annotation table
        enzyme_annot_df = motif_annot.reindex(
            columns=user_inputs["annot_motifs_by"]
        ).sort_values(
            by=user_inputs["annot_motifs_by"]
        )

        # Only keep the motifs which are detected in >=1 genome
        enzyme_annot_df = enzyme_annot_df.reindex(
            index=[
                i for i in enzyme_annot_df.index.values
                if i in plot_value_df.columns.values
            ]
        )

        # Reorder the display data to match
        plot_value_df = plot_value_df.reindex(
            columns=enzyme_annot_df.index.values
        )
        plot_text_df = plot_text_df.reindex(
            columns=enzyme_annot_df.index.values
        )

    # Otherwise, the annotations should be made to match the order of the genomes
    else:
        enzyme_annot_df = motif_annot.reindex(
            columns=user_inputs["annot_motifs_by"] if len(user_inputs["annot_motifs_by"]) > 0 else ['none'],
            index=plot_value_df.columns.values
        )

    # If the user wants to sort the genomes by an annotation
    if user_inputs["sort_genomes_by"] == "Genome Annotations":

        assert len(user_inputs["annot_genomes_by"]) > 0, "Must specify genome annotations for sorting"

        # Sort the genome annotation table
        genome_annot_df = genome_annot.reindex(
            columns=user_inputs["annot_genomes_by"]
        ).sort_values(
            by=user_inputs["annot_genomes_by"]
        )

        # Only keep the genomes which have >=1 motif detected
        genome_annot_df = genome_annot_df.reindex(
            index=[
                i for i in genome_annot_df.index.values
                if i in plot_value_df.index.values
            ]
        )

        # Reorder the display data to match
        plot_value_df = plot_value_df.reindex(
            index=genome_annot_df.index.values
        )
        plot_text_df = plot_text_df.reindex(
            index=genome_annot_df.index.values
        )

    # Otherwise, the annotations should be made to match the order of the genomes
    else:
        genome_annot_df = genome_annot.reindex(
            columns=user_inputs["annot_genomes_by"] if len(user_inputs["annot_genomes_by"]) > 0 else ['none'],
            index=plot_value_df.index.values
        )

    # For the colors, convert all values to numeric and scale to 0-1
    enzyme_marginal_z = enzyme_annot_df.apply(
        convert_text_to_scalar
    )
    genome_marginal_z = genome_annot_df.apply(
        convert_text_to_scalar
    )

    # Set the fraction of the plot used for the marginal annotation
    # depending on the number of those annotations
    enzyme_annot_frac = min(0.5, 0.02 + (0.05 * float(len(user_inputs["annot_motifs_by"]))))
    genome_annot_frac = min(0.5, 0.02 + (0.05 * float(len(user_inputs["annot_genomes_by"]))))

    # If the genomes are being displayed on the horizontal axis
    if user_inputs['genome_axis'] == "Columns":

        # Transpose the DataFrames with genome/motif values
        plot_value_df = plot_value_df.T
        plot_text_df = plot_text_df.T

        # The enzyme marginal annotation will be on the rows
        row_marginal_x = user_inputs["annot_motifs_by"]
        row_marginal_y = plot_value_df.index.values
        row_marginal_z = enzyme_marginal_z.values
        row_marginal_text = enzyme_annot_df.values

        # The genome marginal annotation will be on the columns
        col_marginal_y = user_inputs["annot_genomes_by"]
        col_marginal_x = plot_value_df.columns.values
        col_marginal_z = genome_marginal_z.T.values
        col_marginal_text = genome_annot_df.T.values

        # Add a third row for the barplot
        nrows=3
        ncols=2

        # The size of the marginal plots is driven by the number of annotations
        column_widths = [enzyme_annot_frac, 1 - enzyme_annot_frac]
        row_heights = [genome_annot_frac, 1 - genome_annot_frac, 0.1]

        # Compute the data for the marginal barplot
        bar_x = plot_value_df.columns.values
        bar_y = (plot_value_df > 0).sum(axis=0)
        bar_text = list(map(lambda i: f"{i[0]:,} motifs detected in {i[1]}", zip(bar_y, bar_x)))
        bar_orientation = "v"

    # Otherwise
    else:

        # The genomes must be displayed on the vertical axis
        assert user_inputs['genome_axis'] == "Rows"

        # The genome/motif data does not need to be transposed

        # The enzyme marginal annotation will be on the columns
        col_marginal_x = plot_value_df.columns.values
        col_marginal_y = user_inputs["annot_motifs_by"]
        col_marginal_z = enzyme_marginal_z.T.values
        col_marginal_text = enzyme_annot_df.T.values

        # The genome marginal annotation will be on the rows
        row_marginal_x = user_inputs["annot_genomes_by"]
        row_marginal_y = plot_value_df.index.values
        row_marginal_z = genome_marginal_z.values
        row_marginal_text = genome_annot_df.values

        # Add a third column for the barplot
        nrows=2
        ncols=3

        # The size of the marginal plots is driven by the number of annotations
        row_heights = [enzyme_annot_frac, 1 - enzyme_annot_frac]
        column_widths = [genome_annot_frac, 1 - genome_annot_frac, 0.1]

        # Compute the data for the marginal barplot
        bar_x = (plot_value_df > 0).sum(axis=1)
        bar_y = plot_value_df.index.values
        bar_text = list(map(lambda i: f"{i[0]:,} motifs detected in {i[1]}", zip(bar_x, bar_y)))
        bar_orientation = "h"

    # Set up the figure
    fig = make_subplots(
        rows = nrows,
        cols = ncols,
        vertical_spacing=0.01,
        horizontal_spacing=0.01,
        start_cell="bottom-left",
        column_widths=column_widths,
        row_heights=row_heights,
        shared_xaxes=True,
        shared_yaxes=True,
    )

    # Add the heatmap to the plot
    fig.append_trace(
        go.Heatmap(
            x=plot_value_df.columns.values,
            y=plot_value_df.index.values,
            z=plot_value_df.values,
            colorscale=user_inputs["heatmap_cpal"],
            text=plot_text_df.values,
            hoverinfo="text",
            colorbar_title="Percent<br>Detection<br>of Motif"
        ),
        row=2,
        col=2
    )

    # Add the marginal annotation on the rows
    fig.append_trace(
        go.Heatmap(
            x=row_marginal_x,
            y=row_marginal_y,
            z=row_marginal_z,
            colorscale=user_inputs["annot_cpal"],
            text=row_marginal_text,
            hoverinfo="text",
            showscale=False,
        ),
        row=2,
        col=1
    )

    # Add the marginal annotation on the columns
    fig.append_trace(
        go.Heatmap(
            x=col_marginal_x,
            y=col_marginal_y,
            z=col_marginal_z,
            colorscale=user_inputs["annot_cpal"],
            text=col_marginal_text,
            hoverinfo="text",
            showscale=False
        ),
        row=1,
        col=2
    )

    # Add the barplot with the number of motifs per genome
    fig.append_trace(
        go.Bar(
            x=bar_x,
            y=bar_y,
            text=bar_text,
            hoverinfo="text",
            orientation=bar_orientation,
            marker_color="blue",
        ),
        row=nrows,
        col=ncols
    )

    # Set up the size of the figure
    fig.update_layout(
        height=user_inputs['figure_height'],
        width=user_inputs['figure_width'],
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)'
    )

    # Return the figure
    return fig


def convert_text_to_scalar(c):
    """Take a column with text or strings and convert to values ranging 0-1."""

    # Get the sorted list of values
    unique_values = c.dropna().drop_duplicates().sort_values()
    
    # Assign each value to an index position
    value_map = pd.Series(dict(zip(unique_values, range(len(unique_values)))))

    # Set the maximum value at 1
    value_map = value_map / value_map.max()

    # Map NaN to 0
    return c.apply(value_map.get)


####################
# START THE SCRIPT #
####################

# FORMAT THE REBASE DATA
value_df, text_df, motif_annot, genome_annot = format_rebase_df(os.getcwd())

# TITLE OF THE APP
st.title('Epi-Pangenome Map')

# SIDEBAR MENUS
user_inputs = dict(
    detail_genome=st.sidebar.selectbox(
        "Show Details For",
        list(value_df.index.values) + ["None"]
    ),
    hidden_genomes=st.sidebar.multiselect(
        "Hide Genomes",
        list(value_df.index.values)
    ),
    hidden_motifs=st.sidebar.multiselect(
        "Hide Motifs",
        list(value_df.columns.values)
    ),
    genome_axis=st.sidebar.radio(
        "Show Genomes On",
        ["Rows", "Columns"],
        index=0
    ),
    sort_genomes_by=st.sidebar.selectbox(
        "Sort Genomes By",
        [
            "Motif Presence/Absence",
            "Genome Annotations"
        ]
    ),
    annot_genomes_by=st.sidebar.multiselect(
        "Annotate Genomes By",
        genome_annot.columns.values,
        []
    ),
    sort_motifs_by=st.sidebar.selectbox(
        "Sort Motifs By",
        [
            "Genome Presence/Absence",
            "Motif Annotations"
        ]
    ),
    annot_motifs_by=st.sidebar.multiselect(
        "Annotate Motifs By",
        motif_annot.columns.values,
        ["type_label"]
    ),
    heatmap_cpal=st.sidebar.selectbox(
        "Heatmap Color Palette",
        px.colors.named_colorscales(),
        index=px.colors.named_colorscales().index("blues")
    ),
    annot_cpal=st.sidebar.selectbox(
        "Annotation Color Palette",
        px.colors.named_colorscales(),
        index=px.colors.named_colorscales().index("bluered")
    ),
    figure_height=st.sidebar.slider(
        "Figure Height",
        min_value=100,
        max_value=1200,
        step=1,
        value=600
    ),
    figure_width=st.sidebar.slider(
        "Figure Width",
        min_value=100,
        max_value=1200,
        step=1,
        value=600
    ),
    min_fraction=st.sidebar.slider(
        "Minimum Percentage per Motif",
        min_value=0,
        max_value=100,
        value=25,
        step=1,
        help="Only show motifs which are present in a single genome at this minimum threshold"
    ),
    min_prevalence=st.sidebar.slider(
        "Minimum Number of Genomes per Motif",
        min_value=1,
        max_value=value_df.shape[0],
        value=1,
        step=1,
        help="Only show motifs which are present in at least this many genomes"
    )
)

# MASK ANY SELECTED ROWS/COLUMNS
plot_value_df = value_df.drop(
    columns=user_inputs['hidden_motifs'],
    index=user_inputs['hidden_genomes']
)

# MASK ANY MOTIFS WHICH DO NOT REACH THE MINIMUM THRESHOLD
plot_value_df = plot_value_df.loc[
    :,
    plot_value_df.max() >= user_inputs['min_fraction']
]

# MASK ANY MOTIFS WHICH ARE NOT FOUND IN THE SPECIFIED NUMBER OF GENOMES
if user_inputs['min_prevalence'] > 1:

    # Remove the motifs which do not meet the threshold
    plot_value_df = plot_value_df.reindex(
        columns=plot_value_df.columns.values[
            (plot_value_df > 0).sum() >= user_inputs['min_prevalence']
        ]
    )

    # Remove the genomes which do not contain any motifs
    plot_value_df = plot_value_df.reindex(
        index=plot_value_df.index.values[
            (plot_value_df > 0).sum(axis=1) > 0
        ]
    )


# SORT THE ROWS/COLUMNS
plot_value_df = sort_table(plot_value_df)

# REALIGN TEXT TABLE TO MATCH VALUES
plot_text_df = text_df.reindex(
    index=plot_value_df.index.values,
    columns=plot_value_df.columns.values
)

# MAKE THE MOTIF SUMMARY PLOT
logging.info(f"Plotting {plot_value_df.shape[0]:,} genomes and {plot_value_df.shape[1]:,} enzymes")
st.write(
    format_display(
        plot_value_df,
        plot_text_df,
        motif_annot,
        genome_annot,
        user_inputs
    )
)

# FUNCTION TO FORMAT A DATAFRAME FOR PLOTTING WITH PLOTLY BARPOLAR
def format_bar_df(gbk, fasta, rebase):

    # Map which annotations go into the same track
    annot_type_dict = dict(
        gene="gene",
        CDS="CDS",
        tRNA="RNA",
        rRNA="RNA",
        ncRNA="RNA",
        tmRNA="RNA"
    )

    # The final set of output will have the same set of columns
    bar_df = pd.DataFrame(
        columns=[
            "start",
            "width",
            "track",
            "color",
            "hover_name",
        ]
    )

    # If we have genbank annotations
    if gbk is not None:

        # Add a column which lets us separate annotation types
        gbk = gbk.assign(
            track = gbk["type"].apply(
                annot_type_dict.get
            )
        )

        # Get the length of each contig
        contig_len = gbk.groupby(
            "record_id"
        )["end"].max()

        # Assign global coordinates using an offset
        contig_offset = contig_len.cumsum() - contig_len
        gbk = gbk.assign(
            start = gbk.apply(
                lambda r: min(r['start'], r['end']) + contig_offset[r['record_id']],
                axis=1
            ),
            width = gbk.apply(
                lambda r: abs(r['start'] - r['end']),
                axis=1
            )
        )

        # Remove any records which span entire contigs
        gbk = gbk.loc[
            gbk.apply(
                lambda r: r['width'] < contig_len[r['record_id']],
                axis=1
            )
        ]

        # Add the annotations to the data used for the plot
        bar_df = pd.concat([
            bar_df,
            gbk.assign(
                color=1.,
                # Hover name
                hover_name=lambda d: d.apply(
                    format_gbk_hover_name,
                    axis=1
                )
            ).reindex(
                columns=bar_df.columns.values
            ).dropna()
        ])

    # If we have motif information and the genome sequence
    if rebase is not None and fasta is not None:

        # Iterate over each motif
        for _, r in rebase.iterrows():

            # Only consider rows which have recognition sequences
            if pd.isnull(r.get("rec_seq")):
                continue

            # If this enzyme is in the list of enzymes to skip
            if r["rec_seq"] in user_inputs["hidden_motifs"]:

                # Skip it
                continue
            
            # Add a track for this particular motif
            bar_df = pd.concat([
                bar_df,
                format_motif_genome_track(
                    r,
                    fasta
                )
            ])

    # If we have the genome sequence
    if fasta is not None:

        # Add the GC content
        # Add a track for this particular motif
        bar_df = pd.concat([
            bar_df,
            format_gc_track(fasta)
        ])


    return bar_df


def expand_reverse_complement(nucl_str):
    # Get the reverse complement
    rc_str = str(Seq(nucl_str).reverse_complement())
    
    # If the sequence is a palendrome
    if rc_str == nucl_str:

        # Only check for it once
        return [(nucl_str, "")]

    # If the sequence is not a palendrom
    else:

        # Check for both the forward and the reverse
        return [
            (nucl_str, " (+)"),
            (rc_str, " (-)")
        ]

def format_motif_genome_track(r, fasta, window_size=1000):
    """Format a summary of the density of a particular motif along the genome."""

    logging.info(f"Calculating frequency of {r.rec_seq}")

    # Populate a list of dicts
    output = []

    # Keep track of the offset while iterating over contigs
    offset = 0

    # Iterate over each contig
    for contig_seq in fasta.values():

        # Add a record for each position where the motif matches
        output.extend([
            # Record the position of the match
            {
                "pos": match.start() + offset,
                "strand": strand
            }
            # Iterate over each of the possible sequences which match this motif
            for nucl_str in expand_ambiguous_nucleotides(r.rec_seq)
            # Check both the forward and reverse strand
            for query_str, strand in expand_reverse_complement(nucl_str)
            # Iterate over each position where the sequence matches
            for match in re.finditer(query_str, contig_seq.upper())
        ])

        # Add to the offset for the next contig
        offset += len(contig_seq)

    # Make a DataFrame
    output = pd.DataFrame(output)

    # Get the index of each window
    output = output.assign(
        window = output.pos.apply(lambda v: int(v / window_size))
    )

    # Get the total number of windows
    n_windows = int(offset / window_size)

    # Calculate the proportion of matches per window
    output = pd.concat([
        pd.DataFrame(dict(
            n=strand_df.groupby("window").apply(len),
        )).reindex(
            index=list(range(n_windows))
        ).fillna(
            0
        ).assign(
            track=f"{r.rec_seq}{strand}",
            strand=strand
        )
        for strand, strand_df in output.groupby("strand")
    ]).reset_index()

    # Format the output for the barplot
    output = output.assign(
        # Transform the `window` index to a `start` and `width`
        start=output.window * window_size,
        width=window_size,
        # The `hover_name` will show the number of matches inside the window
        hover_name=output.apply(
            lambda r: f"{int(r.n):,} {'matches' if r.n > 1 else 'match'} in a {window_size:,}bp window",
            axis=1
        ),
        # Scale the color to the maximum for this motif
        color=output.n / output.n.max()
    ).query(
        "n > 0"
    )

    return output.reindex(
        columns=[
            "start",
            "width",
            "hover_name",
            "color",
            "track"
        ]
    )


@st.cache
def format_gc_track(fasta, window_size=1000):
    """Format a summary of GC content along the genome."""

    logging.info(f"Calculating GC content")

    # Populate a list of dicts
    output = []

    # Keep track of the offset while iterating over contigs
    offset = 0

    # Iterate over each contig
    for contig_seq in fasta.values():

        # Add a record for each position where the motif matches
        output.extend([
            # Record the position of the match
            {
                "pos": match.start() + offset,
            }
            # Iterate over each position where the sequence matches
            for match in re.finditer("[GC]", contig_seq.upper())
        ])

        # Add to the offset for the next contig
        offset += len(contig_seq)

    # Make a DataFrame
    output = pd.DataFrame(output)

    # Get the index of each window
    output = output.assign(
        window = output.pos.apply(lambda v: int(v / window_size))
    )

    # Get the total number of windows
    n_windows = int(offset / window_size)

    # Calculate the proportion of matches per window
    output = pd.concat([
        pd.DataFrame(dict(
            n=output.groupby("window").apply(len),
        )).reindex(
            index=list(range(n_windows))
        ).fillna(
            0
        )
    ]).reset_index()

    # Format the output for the barplot
    output = output.assign(
        # Transform the `window` index to a `start` and `width`
        start=output.window * window_size,
        width=window_size,
        # The `hover_name` will show the GC content inside the window
        hover_name=output.apply(
            lambda r: f"{round(100. * r.n / window_size, 1)}% GC",
            axis=1
        ),
        # Scale the color to the maximum for this motif
        color=output.n / output.n.max(),
        track="GC Content"
    )

    return output.reindex(
        columns=[
            "start",
            "width",
            "hover_name",
            "color",
            "track"
        ]
    )

def expand_ambiguous_nucleotides(rec_seq):
    """For a sequence with ambiguous nucleotides, find all the possible exact matches."""

    # Make sure the sequence is uppercase
    rec_seq = rec_seq.upper()

    # Make a list of lists with each of the possibilities at each position
    all_possible = [
        list("." if i == "N" else ambiguous_dna_values[i])
        for i in rec_seq
    ]

    # Make a list of index positions
    for i in product(*all_possible):
        yield ''.join(i)


def format_gbk_hover_name(r):
    """Format the hover text for each annotation."""

    # If there is no locus_tag
    if pd.isnull(r.__dict__.get("locus_tag")):

        # Show the contig and position
        hover_name = f"{r.record_id}: {r.start:,} - {r.end:,}"

    # If there IS a locus tag
    else:

        # Show the locus, contig, and position
        hover_name = f"{r.locus_tag}<br>{r.record_id}: {r.start:,} - {r.end:,}"

    # Additional annotations
    for k in ['product', 'product_id', 'db_xref']:
        if pd.notnull(r.get(k)):
            hover_name = f"{hover_name}<br>{r[k]}"

    return hover_name


# FUNCTION TO RENDER THE GENOME ANNOTATION PLOT
def format_genome_map(genome_name, dirpath):

    # READ INPUT DATA FROM WORKING DIRECTORY
    data, motif_annot, genome_annot = read_data(dirpath)

    # Format a DataFrame which displays the annotations in concentric circles
    bar_df = format_bar_df(
        data[genome_name].get('gbk'),
        data[genome_name].get('fasta'),
        data[genome_name].get('rebase'),
    )

    fig = go.Figure(
        go.Bar(
            base=bar_df["start"],
            x=bar_df["width"],
            y=bar_df["track"],
            marker_color=bar_df["color"],
            marker_line_width=0,
            text=bar_df["hover_name"],
            hoverinfo="text",
            orientation="h",
            marker_colorscale="blues",
        )   
    )

    fig.update_layout(
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        width=700,
        height=700,
        title=genome_name,
    )

    return fig

# MAKE THE GENOME ANNOTATION PLOT
if user_inputs['detail_genome'] != "None":
    logging.info(f"Plotting {user_inputs['detail_genome']}")
    st.write(
        format_genome_map(
            user_inputs['detail_genome'],
            os.getcwd()
        )
    )
else:
    logging.info("User has opted to omit the genome detail display")

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
from plotly.subplots import make_subplots
from scipy import cluster
import streamlit as st
from zipfile import ZipFile


###############################
# DEFINE ALL HELPER FUNCTIONS #
###############################
def read_data(
    dir_path
):
    """Read in data contained within a directory."""

    # Populate a dictionary with the output
    output = defaultdict(dict)

    # By default, there are no additional annotations for the motifs
    addl_motif_annot = None

    # Walk the directory hierarchy
    for dirpath, dirnames, filenames in os.walk(dir_path):

        # Iterate over each file
        for filename in filenames:

            # If this is a CSV, then it may be the motif annotation table
            if filename.endswith('.csv'):

                # Try to read it in
                try:
                    addl_motif_annot = pd.read_csv(
                        os.path.join(dirpath, filename),
                        index_col=0
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
            ]
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
    ).set_index("rec_seq")

    # Now add in any additional annotations provided by the user
    if addl_motif_annot is not None:
        
        # Iterate over the columns of data provided by the user
        for col_name, col_values in addl_motif_annot.items():
            
            # Assign the values provided by the user to the rebase data
            motif_annot = motif_annot.assign(
                **{
                    col_name: col_values
                }
            )

    # Finally, reformat the index used to label the motif
    # so that it includes both the enzyme type and subtype
    # as indicated by the rebase data
    motif_annot = motif_annot.rename(
        index={
            r['rec_seq']: format_enzyme_name(r)
            for _, r in motif_annot.reset_index().iterrows()
        }
    )

    # Return the dict of all data, along with the additional motif annotations, if any
    return output, motif_annot


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
                data[-1][k] = v

    return pd.DataFrame(data)

def parse_rebase(lines):
    """Parse the epigenetic data from REBASE."""

    # Populate a list with dicts, each entry being an enzyme
    enzymes = list()
    enzyme = dict()

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

                # Add it to the list
                enzymes.append(enzyme)

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
    data, motif_annot = read_data(dirpath)

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
        enzyme_name = df.apply(
            format_enzyme_name,
            axis=1
        ),
        text = df.apply(
            format_rebase_text,
            axis=1
        )
    )

    # Transform to a wide DataFrame containing floats for `percent_detected`
    percent_detected = df.pivot(
        index="organism",
        columns="enzyme_name",
        values="percent_detected"
    ).fillna(
        0
    ).applymap(
        float
    )

    # Also make a wide table with the complete set of REBASE outputs
    text_df = df.pivot(
        index="organism",
        columns="enzyme_name",
        values="text"
    ).fillna(
        ""
    ).reindex(
        index=percent_detected.index.values,
        columns=percent_detected.columns.values,
    )

    return percent_detected, text_df, motif_annot


def format_enzyme_name(r):
    """Format the label of each motif from the REBASE TXT file."""

    # All enzymes must have a rec_seq
    assert not pd.isnull(r.get("rec_seq")), r

    # Format the base name
    enzyme_name = r.rec_seq

    # If the enzyme type is present
    if not pd.isnull(r.get("enz_type")):

        # Add it to the label
        enzyme_name = f"{enzyme_name} - Type {int(r.enz_type)}"

        # If the subtype is present
        if not pd.isnull(r.get("sub_type")):

            # Add it to the label
            enzyme_name = f"{enzyme_name}{r.sub_type}"

    return enzyme_name


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

def format_enzyme_type(enzyme_name_list):
    """Format the data to display a marginal annotation of enzyme type."""
    
    enzyme_type_list = [
        f"Type {enzyme_name.rsplit(' - Type ', 1)[1]}" if " - Type " in enzyme_name else None
        for enzyme_name in enzyme_name_list
    ]

    # Convert to numeric
    enzyme_type_mapping = {
        enzyme_type: enzyme_type if enzyme_type is None else i
        for i, enzyme_type in enumerate(
            pd.Series(
                enzyme_type_list
            ).drop_duplicates(
            ).sort_values(
            )
        )
    }
    enzyme_ix_list = list(map(enzyme_type_mapping.get, enzyme_type_list))

    return pd.DataFrame(
        dict(
            enzyme_type=enzyme_type_list,
            ix=enzyme_ix_list,
        ),
        index=enzyme_name_list
    )

def format_display(plot_value_df, plot_text_df, motif_annot, user_inputs):
    """Format the display."""

    # Format the enzyme type key, which includes both the nicely-formatted
    # text describing the enzyme type and subtype, as well as a numeric
    # index which helps provide a color assignment in the plot
    enzyme_annot_df = format_enzyme_type(
        plot_value_df.columns.values
    )

    # Combine that table with any of the additional motif annotations available
    for col_name, col_values in motif_annot.items():
        enzyme_annot_df = enzyme_annot_df.assign(
            **{
                col_name: col_values
            }
        )

    # If the user wants to sort the enzymes by an annotation
    if len(user_inputs["sort_enzymes_by"]) > 0:

        # Sort the annotation table
        enzyme_annot_df.sort_values(
            by=user_inputs["sort_enzymes_by"],
            inplace=True
        )

        # Reorder the display data to match
        plot_value_df = plot_value_df.reindex(
            columns=enzyme_annot_df.index.values
        )
        plot_text_df = plot_text_df.reindex(
            columns=enzyme_annot_df.index.values
        )

    # Set up a list of columns to annotate enzymes by
    marginal_z_columns = ["ix"] + user_inputs["annot_enzymes_by"]
    marginal_text_columns = ["enzyme_type"] + user_inputs["annot_enzymes_by"]
    marginal_labels = ["Enzyme Type"] + user_inputs["annot_enzymes_by"]

    # Set up the annotations which will be displayed in the margin
    marginal_text = enzyme_annot_df.reindex(columns=marginal_text_columns)

    # For the colors, convert all values to numeric and scale to 0-1
    marginal_z = enzyme_annot_df.reindex(
        columns=marginal_z_columns
    ).applymap(
        lambda v: pd.to_numeric(v, errors="coerce")
    ).apply(
        lambda c: c - c.min() if pd.notnull(c.min()) else c
    ).apply(
        lambda c: c / c.max() if pd.notnull(c.max()) else c
    )

    # Set the fraction of the plot used for the marginal annotation
    # depending on the number of those annotations
    annot_frac = 0.02 + (0.05 * float(len(marginal_labels)))

    # If the genomes are being displayed on the horizontal axis
    if user_inputs['genome_axis'] == "Columns":

        # The subplots will be laid out side-by-side
        nrows=1
        ncols=2
        domains=dict(
            xaxis_domain=[annot_frac + 0.01, 1.0],
            xaxis2_domain=[0, annot_frac - 0.01],
            yaxis_domain=[0, 1.0],
            yaxis_anchor="x2"
        )

        # Transpose the DataFrames
        plot_value_df = plot_value_df.T
        plot_text_df = plot_text_df.T

        # The marginal annotation will be vertical
        marginal_x = marginal_labels
        marginal_y = plot_value_df.index.values
        marginal_z = marginal_z.values
        marginal_text = marginal_text.values

    # Otherwise
    else:

        # The genomes must be displayed on the vertical axis
        assert user_inputs['genome_axis'] == "Rows"

        # The subplots will be laid out top-and-bottom
        nrows=2
        ncols=1
        domains=dict(
            yaxis_domain=[annot_frac + 0.01, 1.0],
            yaxis2_domain=[0, annot_frac - 0.01],
            xaxis_domain=[0, 1.0],
            xaxis_anchor="y2"
        )

        # The data does not need to be transposed

        # The marginal annotation will be horizontal
        marginal_x = plot_value_df.columns.values
        marginal_y = marginal_labels
        marginal_z = marginal_z.T.values
        marginal_text = marginal_text.T.values

    # Set up the figure
    fig = make_subplots(
        rows = nrows,
        cols = ncols
    )

    # Add the heatmap to the plot
    fig.append_trace(
        go.Heatmap(
            x=plot_value_df.columns.values,
            y=plot_value_df.index.values,
            z=plot_value_df.values,
            colorscale=user_inputs["heatmap_cpal"],
            text=plot_text_df.values,
            hoverinfo="text"
        ),
        row=1,
        col=1
    )

    # Add the marginal annotation of enzyme type
    fig.append_trace(
        go.Heatmap(
            x=marginal_x,
            y=marginal_y,
            z=marginal_z,
            colorscale=user_inputs["annot_cpal"],
            text=marginal_text,
            hoverinfo="text",
            showscale=False,
        ),
        row=nrows,
        col=ncols
    )

    # Set up the size of the figure
    fig.update_layout(
        height=user_inputs['figure_height'],
        width=user_inputs['figure_width'],
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        **domains
    )

    # Return the figure
    return fig


####################
# START THE SCRIPT #
####################

# FORMAT THE REBASE DATA
value_df, text_df, motif_annot = format_rebase_df(os.getcwd())

# TITLE OF THE APP
st.title('Epi-Pangenome Map')

# SIDEBAR MENUS
user_inputs = dict(
    hidden_genomes=st.sidebar.multiselect(
        "Hide Genomes",
        list(value_df.index.values)
    ),
    hidden_enzymes=st.sidebar.multiselect(
        "Hide Enzymes",
        list(value_df.columns.values)
    ),
    genome_axis=st.sidebar.radio(
        "Show Genomes On",
        ["Rows", "Columns"],
        index=0
    ),
    sort_enzymes_by=st.sidebar.multiselect(
        "Sort Enzymes By",
        motif_annot.columns.values
    ),
    annot_enzymes_by=st.sidebar.multiselect(
        "Annotate Enzymes By",
        motif_annot.columns.values
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
    columns=user_inputs['hidden_enzymes'],
    index=user_inputs['hidden_genomes'],
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

# MAKE THE PLOT
logging.info(f"Plotting {plot_value_df.shape[0]:,} genomes and {plot_value_df.shape[1]:,} enzymes")
st.write(
    format_display(
        plot_value_df,
        plot_text_df,
        motif_annot,
        user_inputs
    )
)

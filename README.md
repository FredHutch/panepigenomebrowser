# panepigenomebrowser
Browser of epigenetic data from multiple bacterial isolates

## Install Prerequisites

To install prerequisites, first set up a python3 virtual environment:

```#!/bin/bash
# Create the virtual environment
python3 -m venv venv

# Activate the virtual environment
source venv/bin/activate

# Install the prerequisites for running this app
# Replacing $REPO_FOLDER with the folder containing this repository
python3 -m pip install -r $REPO_FOLDER/requirements.txt
```

## Launch the app

To render data from a set of files in the local working directory,
run the following command:

```#!/bin/bash

# Replace $REPO_FOLDER with the folder containing this repository
streamlit run $REPO_FOLDER/app.py
```

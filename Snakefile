""" Snake pipeline for running cellranger with CITE-seq data """

# This has been borrowed and modified from pipelines written by Ryan Sheridan

# Configure shell for all rules 
shell.executable("/bin/bash")
shell.prefix("set -o nounset -o pipefail -o errexit -x; ")
import subprocess
import glob
import os 
import re
from collections import defaultdict



# Parameters from config.yaml
RAW_DATA        = config["RAW_DATA"]
SAMPLES         = config["SAMPLES"]
RNA_SAMPLES     = config["RNA_SAMPLES"]
ADT_SAMPLES     = config["ADT_SAMPLES"]
VDJ_T_SAMPLES   = config["VDJ_T_SAMPLES"]
VDJ_B_SAMPLES   = config["VDJ_B_SAMPLES"]
RESULTS         = config["RESULTS"]
GENOME          = config["GENOME"]
ADT_REF         = config["ADT_REF"]
VDJ_REF         = config["VDJ_REF"]
MAX_JOBS        = config["MAX_JOBS"]
LSF_TEMPLATE    = config["LSF_TEMPLATE"]
AGGR_GROUP      = config["AGGR_SAMPLES"]
CHEMISTRY       = config["CHEMISTRY"]
VELOCYTO_GROUP  = config["VELOCYTO_GROUP"]
SCRIPT_PATH     = config["SCRIPT_PATH"]
SCRIPTS_RUN     = config["SCRIPTS_RUN"]
SAMPLE_SCRIPTS  = config["SAMPLE_SCRIPTS"]
MERGE_SCRIPTS   = config["MERGE_SCRIPTS"]
SAMPLE_METADATA = config["SAMPLE_METADATA"]
SAMPLE_INFO     = config["SAMPLE_INFO"]

# Function to check paths for input files/directories
def _check_path(path):
    if os.path.exists(path):
        return os.path.abspath(path)
    else:
        sys.exit("ERROR: " + path + " does not exist.")



# Set sample/group names
RNA_SAMPLES     = [x.strip() for x in RNA_SAMPLES]
FASTQ_INFO      = "_S[0-9]+_L[0-9]+_R[12]_[0-9]+\.fastq\.gz"
FASTQ_SHORT     = "_S[0-9]+_R[12]_[0-9]+\.fastq\.gz"

SAMPLE_DICT_RNA = {SAMPLES[i]: RNA_SAMPLES[i] for i in range(len(SAMPLES))}

if ADT_SAMPLES:
    ADT_SAMPLES = [re.sub(", ", ",", x.strip()) for x in ADT_SAMPLES]
    #SAMPLES     = [x + "-" + re.sub(",", "_", y) for x, y in zip(RNA_SAMPLES, ADT_SAMPLES)]
    SAMPLE_DICT_ADT = {SAMPLES[i]: ADT_SAMPLES[i] for i in range(len(SAMPLES))}
    ADT_REF     = _check_path(ADT_REF)
else:
    SAMPLE_DICT_ADT = {SAMPLES[i]: "" for i in range(len(SAMPLES))}

if VDJ_T_SAMPLES:
    VDJ_T_SAMPLES = [x.strip() for x in VDJ_T_SAMPLES]
    SAMPLE_DICT_VDJ_T = {SAMPLES[i]: VDJ_T_SAMPLES[i] for i in range(len(SAMPLES))}
    VDJ_REF     = _check_path(VDJ_REF)
else:
    SAMPLE_DICT_VDJ_T = {SAMPLES[i]: "" for i in range(len(SAMPLES))}


if VDJ_B_SAMPLES:
    VDJ_B_SAMPLES = [x.strip() for x in VDJ_B_SAMPLES]
    SAMPLE_DICT_VDJ_B = {SAMPLES[i]: VDJ_B_SAMPLES[i] for i in range(len(SAMPLES))}
    VDJ_REF     = _check_path(VDJ_REF)
else:
    SAMPLE_DICT_VDJ_B = {SAMPLES[i]: "" for i in range(len(SAMPLES))}

# Make a chemistry dictionary
for i in range(len(SAMPLES)):
    if(CHEMISTRY is not None):
        if not CHEMISTRY.get(SAMPLES[i]):
            CHEMISTRY[SAMPLES[i]] = "auto"
    else:
        CHEMISTRY = {}
        CHEMISTRY[SAMPLES[i]] = "auto"

# Check/set directory/file paths
if len(RAW_DATA) == len(SAMPLES):
    RAW_DATA_DICT = {SAMPLES[i]: RAW_DATA[i] for i in range(len(SAMPLES))}
elif len(RAW_DATA) == 1:
    RAW_DATA_DICT = {SAMPLES[i]: RAW_DATA[0] for i in range(len(SAMPLES))}
else:
    sys.exit("RAW_DATA must either be the same length of samples or length 1")

if not os.path.exists(RESULTS):
    os.makedirs(RESULTS)
RESULTS = _check_path(RESULTS)
GENOME  = _check_path(GENOME)

FASTQ_DIR = RESULTS + "/fastqs"
if not os.path.exists(FASTQ_DIR):
    os.makedirs(FASTQ_DIR)

if LSF_TEMPLATE:
    LSF_TEMPLATE = _check_path(LSF_TEMPLATE)
else:
    LSF_TEMPLATE = "lsf"

if not AGGR_GROUP:
    AGGR_GROUP = "none"

if not VELOCYTO_GROUP:
    VELOCYTO_GROUP = "none"

if SCRIPTS_RUN:
    # Get R script info
    # All script names is the keys of scripts run
    script_name_dict = {i:"all_scripts" for i in SCRIPTS_RUN.keys()}


    # Create script dict where all scripts is all the scripts_run scripts
    script_dict = {"all_scripts": SCRIPTS_RUN}
else:
    script_dict = {}
    script_name_dict = {}

if SAMPLE_SCRIPTS:
    if SAMPLE_SCRIPTS["samples"] == "all":
        SAMPLE_SCRIPTS["samples"] = SAMPLES
    for sample_name in SAMPLE_SCRIPTS["samples"]:
        individual_dict = OrderedDict()
        for script in SAMPLE_SCRIPTS["scripts_run"]:
            # Get key and value information
            dict_key = sample_name + "__" + script
            dict_key = re.sub('\\.R', '', dict_key)
            dict_value = os.path.join("indiviual_analysis", script)
            # Add to an individual dict
            individual_dict[dict_key] = dict_value
            # Add the script name to all names
            script_name_dict[dict_key] = sample_name
        # Add new dict to the script dict under the sample name as the key
        script_dict[sample_name] = individual_dict

if MERGE_SCRIPTS:
    if MERGE_SCRIPTS["samples"] == "all":
        MERGE_SCRIPTS["samples"] = SAMPLES
    merged_dict = OrderedDict()
    for script in MERGE_SCRIPTS["scripts_run"]:
        dict_key = "merged__" + script
        dict_key = re.sub('\\.R', '', dict_key)
        dict_value = os.path.join("integrated_analysis", script)
        merged_dict[dict_key] = dict_value
        script_name_dict[dict_key] = "merged"
    script_dict["merged"] = merged_dict

# Check that all names were entered correctly
for name in SAMPLES:
    name_1 = re.search("_[A-Z]+", name).group()
    name_2 = re.search("_[A-Z]+", SAMPLE_DICT_ADT[name]).group()
    name_3 = re.search("_[A-Z]+", SAMPLE_DICT_RNA[name]).group()
    name_4 = re.search("_[A-Z]+", SAMPLE_DICT_VDJ_B[name]).group()
    if(name_1 != name_2):
        sys.exit(name_1 + " does not equal " + name_2 + " in the ADT list.")
    if(name_1 != name_3):
        sys.exit(name_1 + " does not equal " + name_2 + " in the RNA list.")
    if(name_1 != name_4):
        sys.exit(name_1 + " does not equal " + name_2 + " in the VDJ_B list.")

# Final output files
rule all:
    input:
        expand(
            "{results}/logs/{sample}_csv_done.txt",
            results = RESULTS, sample = SAMPLES
        ),
        expand(
            "{results}/logs/{sample}_count_done.txt",
            results = RESULTS, sample = SAMPLES
            ),
        expand(
            "{results}/logs/{group}_csv_aggr_done.txt",
            results = RESULTS, group = AGGR_GROUP
            ),
        expand(
            "{results}/logs/{group}_cellranger_aggr_done.txt",
            results = RESULTS, group = AGGR_GROUP
            ),
        # Run dropkick to identify cells
        expand(
            "{results}/R_analysis/{sample}/files/dropkick_cells.csv",
            results = RESULTS, sample = SAMPLES
            ),
        expand(
            "{results}/R_analysis/{sample}/files/scar_denoised.csv",
            results = RESULTS, sample = SAMPLES
            ),
        expand(
            "{results}/run_scripts/{sample}_run.txt",
            sample = script_name_dict.keys(), results = RESULTS
            ),
        # Add in immcantation following 
        # https://immcantation.readthedocs.io/en/stable/tutorials/10x_tutorial.html
        expand(
            "{results}/{sample}/outs/immcantation/{sample}_finished.txt",
            results = RESULTS, sample = SAMPLES
            ),
        expand(
            "{results}/R_analysis/merged/define_clones/all_clones_finished.txt",
            results = RESULTS
            ),
        expand(
            "{results}/t1k/{sample}_t1k_finished.txt",
            results = RESULTS, sample = SAMPLES
            ),
        expand(
             "{results}/t1k/combined_results.tsv",
             results = RESULTS
        )
        # Add in germlines
        #expand(
        #    "{results}/R_analysis/merged/define_clones/germlines_finished.txt",
        #    results = RESULTS
        #    )
        # expand(
        #     "{results}/logs/{sample}_velocyto_done.out",
        #     results = RESULTS, sample = SAMPLES
        # ),
        # expand(
        #     "{results}/logs/{group}_velocyto_combined.out",
        #     results = RESULTS, group = VELOCYTO_GROUP
        #     )

include: "src/rules/cellranger_multi.snake"
include: "src/rules/velocyto.snake"
include: "src/rules/run_rscripts.snake"
include: "src/rules/immcantation.snake"
include: "src/rules/run_t1k.snake"

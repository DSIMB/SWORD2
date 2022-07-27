#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import argparse
import logging
import math
import multiprocessing
import os
import random
import re
import shlex
import shutil
import plotly.express as px
import subprocess
import sys
import time
from copy import copy
from functools import partial

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from matplotlib import patches
from prody import *
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry


def check_parsing_pdb(uniprot_code, pdb_code, pdb_chain, pdb_file):
    """
    This function tries to fetch and parse the input PDB submitted,
    either PDB code and chain or a whole PDB file downloaded by the user.

    Args:
        - uniprot_code (str): AlphaFold Uniprot Accession Id
        - pdb_code (str): PDB code to fetch and parse
        - pdb_chain (str): PDB chain to fetch and parse
        - pdb_file (str): The path to the file downloaded by the user


    Returns:
        - prot (ProDy Protein object): if fetched and parsed correctly, else False
    """
    # Custom user file
    if pdb_file:
        logging.info("Try to parse user input PDB file")
        try:
            prot = parsePDB(pdb_file, chain=pdb_chain)
        except Exception as e:
            sys.exit(str(e))
        if prot is None:
            # Try again without specifying the chain. If it works and
            # if the protein has no chain, we set it to "A" by default.
            # It is necessary because the scoring program
            # accepts only PDB files which contain a chain
            prot = parsePDB(pdb_file)
            if prot is not None and prot.getChids()[0] == ' ':
                prot.setChids([pdb_chain for i in range(prot.numAtoms())])
            else:
                sys.exit(f"Atomic data could not be parsed for chain {pdb_chain}. Please check the input PDB file")
        # Clean non-standard aa
        prot = prot.select('protein and not nonstdaa')
        if prot is None:
            sys.exit(f"Error: No atomic data is left after trying to keep the 20 classical residues. Please check your PDB file.")
    # Download the AlphaFold model
    elif uniprot_code:
        logging.info(f"Download AlphaFold Uniprot Accession ID: {uniprot_code}")
        ok, response = download_af_model(uniprot_code)
        if not ok:
            sys.exit(f"Error: {response}. Please try again.")
        else:
            pdb_file = response
        # Redirect useful error of ProDy
        try:  # Try to parse PDB file
            prot = parsePDB(pdb_file, chain=pdb_chain)
        except Exception as e:
            sys.exit(str(e))
        # A parsed PDB by ProDy returns an AtomGroup, else it could be an EMD file, which we don't want...
        if type(prot) is not AtomGroup:
            sys.exit(f"Error: Something went wrong with the AlphaFold Uniprot Accession Id {uniprot_code}. Please make sure chain {pdb_chain} is valid and that it is a valid ID referenced in the AlphaFold database (https://alphafold.ebi.ac.uk/)")
        if prot is None:
            sys.exit(f"Error: Atomic data could not be parsed. Please check the PDB file corresponding to the code {pdb_code}. Also check that it actually contains the chain {pdb_chain}.")
        # Clean non-standard aa
        prot = prot.select('protein and not nonstdaa')
        if prot is None:
            sys.exit(f"Error: No atomic data is left after trying to keep the 20 classical residues. Please check your PDB file.")
        pdb_file = os.path.basename(pdb_file)
    # Fetch and parse a PDB from a given PDB code and chain
    else:
        logging.info(f"Fetch PDB ID: {pdb_code}")
        # Redirect useful error of ProDy
        try:  # Try to parse fetched PDB file
            prot = parsePDB(pdb_code, chain=pdb_chain, folder=RESULTS_DIR, compressed=False)
        except Exception as e:
            sys.exit(str(e))
        # A parsed PDB by ProDy returns an AtomGroup, else it could be an EMD file, which we don't want...
        if type(prot) is not AtomGroup:
            sys.exit(f"Error: No PDB file could be parsed. Please check that the PDB code {pdb_code} exists in the PDB RCSB database (https://www.rcsb.org) with a legacy PDB format file available, and contains the chain {pdb_chain}. Careful, 'A' is different than 'a'. Please note that mmCIF files are not yet supported.")
        if prot is None:
            sys.exit(f"Error: Atomic data could not be parsed. Please check the PDB file corresponding to the code {pdb_code}. Also check that it actually contains the chain {pdb_chain}.")
        # Clean non-standard aa
        prot = prot.select('protein and not nonstdaa')
        if prot is None:
            sys.exit(f"Error: No atomic data is left after trying to keep the 20 classical residues. Please check your PDB file.")
    return prot

def requests_retry_session(retries=3,
                           backoff_factor=0.3,
                           status_forcelist=(500, 502, 504),
                           session=None):
    session = session or requests.Session()
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session

def download_af_model(id):
    """
    Download the Alphafold2 model corresponding to the Uniprot Id given by user
    https://alphafold.ebi.ac.uk/

    Returns:
        - File path (string): Path of the downloaded PDB file
        or
        False if wrong id
        "DOWNLOAD ERROR" if could not download
    """
    content = []
    with open(f"{BASE_DIR}/data/accession_ids.txt", "r") as f:
        for line in f:
            content.append(line.strip().split(","))
    uniprot_to_infos_dict = {uniprot_acc: [first_res, last_res, af_db_id, latest_version] for uniprot_acc, first_res, last_res, af_db_id, latest_version in content}
    try:
        af_id = uniprot_to_infos_dict[id][2]
        version = uniprot_to_infos_dict[id][3]
    except KeyError:
        sys.exit(f"{id} was not found in AlphaFold DB")
    name = f"{af_id}-model_v{version}"
    url = f"https://alphafold.ebi.ac.uk/files/{name}.pdb"
    try:
        response = requests_retry_session().get(url)
    except Exception as x:
        sys.exit(x)
    with open(f"{RESULTS_DIR}/{name}.pdb", "w") as f:
        f.write(response.text)
    return (True, f"{RESULTS_DIR}/{name}.pdb")


def parse_sword(output):
    """
    Parse the SWORD output and return a dictionary.

    Args:
        - output (list of str): The output of SWORD as a list of strings

    Returns:
        - sword_results (dict): Dictionary containing all partitioning assignments
                                made by SWORD and Peeling
    """
    logging.info("Parse SWORD output")
    # Index of the alternative partitionings
    nb_alt = 0
    # Dict containing all the informations given by SWORD2 output
    sword_results = {"DOMAINS": {}, "AMBIGUITY": "n/a"}
    for line in output:
        amb_found = re.search(r"^A-index = (\++)$", line)
        ass_found = re.search(r"\d{1,}\s+\|", line)
        # Found Ambiguity-index line
        if amb_found:
            sword_results["AMBIGUITY"] = amb_found.group(1)
        # Found a domain assignment
        elif ass_found:
            splitted_ass = re.split(r"\s{0,}\|\s{0,}", line)
            sword_results["DOMAINS"][nb_alt] = {}
            sword_results["DOMAINS"][nb_alt]["NB_DOMAINS"] = int(splitted_ass[0].strip())
            sword_results["DOMAINS"][nb_alt]["MIN_SIZE"] = int(splitted_ass[1].strip())
            boundaries = re.split(r"\s", splitted_ass[2].strip())
            sword_results["DOMAINS"][nb_alt]["BOUNDARIES"] = {}
            for i, boundary in enumerate(boundaries):
                mult_boundaries = re.split(r";", boundary)
                sword_results["DOMAINS"][nb_alt]["BOUNDARIES"][i] = []
                for mb in mult_boundaries:
                    start_pu, end_pu = re.split(r"-", mb)
                    sword_results["DOMAINS"][nb_alt]["BOUNDARIES"][i].append((int(start_pu), int(end_pu)))
            sword_results["DOMAINS"][nb_alt]["AVERAGE K"] = float(splitted_ass[3].strip())
            sword_results["DOMAINS"][nb_alt]["QUALITY"] = splitted_ass[4].strip()
            nb_alt += 1
    return sword_results


def get_quality_as_nb_bars(quality):
    """

    """
    nb_bars = 0
    if quality == "*****" or quality == "+++++":
        nb_bars = 5
    if quality == "****" or quality == "++++":
        nb_bars = 4
    elif quality == "***" or quality == "+++":
        nb_bars = 3
    elif quality == "**" or quality == "++":
        nb_bars = 2
    elif quality == "*" or quality == "+":
        nb_bars = 1
    elif quality == "n/a":
        nb_bars = 0
    return nb_bars

def write_partitionings(sword_results, energies):
    """
    Write the partitionnings into HTML formated text file
    that will be loaded into the results web page.

    Args:
        - sword_results (dict): Dictionary containing all partitioning assignments
                                made by SWORD and Peeling
    Returns:
        - None
    """
    logging.info("Write the SWORD results")
    partitioning = os.path.join(RESULTS_DIR, "SWORD2_summary.txt")
    with open(partitioning, "w") as f:
        # AMBIGUITY INDEX
        nb_bars = get_quality_as_nb_bars(sword_results["AMBIGUITY"])
        f.write("Ambiguity index: " + "*" * nb_bars + "\n")
        for nb_alt_part, alt_part in sword_results["DOMAINS"].items():
            f.write(f"""-----------------------\n""")
            # OPTIMAL PARTITIONING
            if nb_alt_part == 0:
                f.write("""Optimal partition\n""")
            # ALTERNATIVE PARTITIONINGS
            else:
                f.write(f"""Alternative partition {nb_alt_part}\n""")
            nb_bars = get_quality_as_nb_bars(alt_part["QUALITY"])
            f.write("Quality: " + "*" * nb_bars + "\n")
            f.write(f"""Nb. domains: {len(alt_part["BOUNDARIES"])}\n""")
            for i, dom in alt_part["BOUNDARIES"].items():
                f.write(f"""Domain:{i+1}       AUL={int((1-(1/(energies[(nb_alt_part, i)][1])**2))*100) if abs(energies[(nb_alt_part, i)][1]) >= 1 else 0:3}% Z-score={round(energies[(nb_alt_part, i)][1], 1)}\n""")
                for j, (start_pu, end_pu) in enumerate(dom):
                    f.write(f"""    PU:{str(start_pu)+"-"+str(end_pu):>7} AUL={int((1-(1/(energies[(nb_alt_part, i, start_pu, end_pu)][1])**2))*100) if abs(energies[(nb_alt_part, i, start_pu, end_pu)][1]) >= 1 else 0:3}% Z-score={round(energies[(nb_alt_part, i, start_pu, end_pu)][1], 1)}\n""")


def define_colors(sword_results):
    """
        Set visually distinct colors for Domains (dark palette) and PUs (ligh palette)

        Return:
            - pus_colors (dict): key=(start_pu, end_pu) --> value=(r, g, b)
            - dom_colors (dict): key=domain_id --> value=(r, g, b)
    """
    pus_colors = {}
    dom_colors = {}
    color_domain_cnt = 0
    color_pu_cnt = 0
    dom_id = None
    # Set of 20 visually distinct colors for PUs (light palette)
    # https://medialab.github.io/iwanthue/
    colors_for_pus = ["#baeae5", "#e1c65b", "#b4bcf7", "#d0e47b", "#f0a8e5", "#6de4ac", "#d8c0e4", "#a5e18d", "#68d1f1", "#f3b175", "#63e3d8", "#ebbaba", "#c3d28c", "#aac5e2", "#e8da92", "#bcdbec", "#e1c298", "#98c7c6", "#abddb4", "#d4d8bb"]

    # Set of 25 visually distinct colors for domains (dark palette)
    # https://medialab.github.io/iwanthue/
    colors_for_domains = ["#27a3b4", "#c08423", "#d83e7c", "#986a35", "#686fdf", "#559a3b", "#763da6", "#8d8d36", "#ce61c7", "#406021", "#a42c88", "#3d956b","#341d79", "#cf3a44", "#3c8cc9", "#cf6430", "#4d4b92", "#7d3119", "#7b81cd", "#cf6c61", "#401d56", "#c95f7a", "#792e65", "#82263a", "#b86da8"]
    for i, part in sword_results["DOMAINS"].items():
        # COLOR ALL PUS OF A SWORD ALTERNATIVE DOMAIN WITH A DIFFERENT COLOR
        for j, dom in part["BOUNDARIES"].items():
            # Consider that a domain is a list of PUs delineation sorted by 1st delimitation of PUs
            dom_id = tuple(sorted(dom, key=lambda x: x[0]))
            if dom_id not in dom_colors:
                # Pick a new color
                hex_val = colors_for_domains[color_domain_cnt].lstrip("#")
                (r, g, b) = tuple(int(hex_val[k:k+2], 16) for k in (0, 2, 4))
                dom_colors[dom_id] = (r, g, b)
                color_domain_cnt += 1
        for j, dom in part["BOUNDARIES"].items():
            dom_id = tuple(sorted(dom, key=lambda x: x[0]))
            r, g, b = dom_colors[dom_id]
        for j, dom in part["BOUNDARIES"].items():
            dom_id = tuple(sorted(dom, key=lambda x: x[0]))
            for start_pu, end_pu in dom:
                if (start_pu, end_pu) not in pus_colors:
                    # Pick a new color
                    hex_val = colors_for_pus[color_pu_cnt].lstrip("#")
                    (r, g, b) = tuple(int(hex_val[k:k+2], 16) for k in (0, 2, 4))
                    pus_colors[(start_pu, end_pu)] = (r, g, b)
                    color_pu_cnt += 1
    return pus_colors, dom_colors


def write_domains_histogram(sword_results, domains_colors):
    """
    Generate histogram of SWORD domains consistency.
    """
    logging.info("Generate histogram of SWORD2 domains consistency")
    histogram = os.path.join(RESULTS_DIR, "domains_histogram.png")
    domains = {}
    for i, part in sword_results["DOMAINS"].items():
        for j, domain in part["BOUNDARIES"].items():
            dom_id = tuple(sorted(domain, key=lambda x: x[0]))
            if dom_id not in domains:
                domains[dom_id] = {}
            if "nb" not in domains[dom_id]:
                domains[dom_id]["nb"] = 1
            else:
                domains[dom_id]["nb"] += 1
    sorted_domains_by_decreasing_count = sorted(domains.items(), key=lambda x: x[1]["nb"], reverse=True)
    x = []
    for dom_part, infos in sorted_domains_by_decreasing_count:
        x.append(str(list(dom_part)).strip("[]"))
    y = []
    for dom_part, infos in sorted_domains_by_decreasing_count:
        y.append(infos["nb"])
    colors = {}
    for dom_part, infos in sorted_domains_by_decreasing_count:
        r, g, b = domains_colors[dom_part]
        colors[dom_part] = f"rgb({r}, {g}, {b})"
    
    # Calling DataFrame constructor after zipping
    # both lists, with columns specified
    df = pd.DataFrame(list(zip(x, y)),
                columns =['SWORD Domains', 'Count'])
    fig = px.bar(df, x="SWORD Domains", y="Count", color=colors, title="Consistency of domains determined by SWORD")
    fig.update_layout(showlegend=False)
    fig.write_image(histogram, scale=4)


def predict_time(prot):
    """ Time predicted in function of len of protein in minutes"""
    return str(int(133*math.exp(2.06*10**-3*len(set(prot.getResnums())))/60))


def get_energy_and_z_score(pdb, res_list=None):
    """
    Calculate pseudo-energy and z-score of a protein or specified residue list
    of the input pdb.
    List is res+chain: 10A,11A,12A,13A
    """
    if res_list:
        cmd_args = f"{BIN_DIR}/mypmfs-master/scoring -i {pdb} -d {BIN_DIR}/mypmfs-master/025_30_100_potential -q {res_list} -z -s 2000"
    else:
        cmd_args = f"{BIN_DIR}/mypmfs-master/scoring -i {pdb} -d {BIN_DIR}/mypmfs-master/025_30_100_potential -z -s 2000"
    cmd_args = shlex.split(cmd_args)
    output = subprocess.run(cmd_args, capture_output=True, check=True)
    output = output.stdout.decode("utf-8")
    output = output.split("\n")
    energy = None
    z_score = None
    for line in output:
        pseudo_e_found = re.search(r"^Pseudo-energy = (.+)$", line)
        z_score_found = re.search(r"^Z-score = (.+)$", line)
        if pseudo_e_found:
            energy = float(pseudo_e_found.group(1))
        if z_score_found:
            z_score = float(z_score_found.group(1))
    return energy, z_score


def multiprocess_get_energy(i, energies, dom_bounds):
    """
    Calculate the energy and Z-score of PUs and Domains.

    Args:
        - i: index of partitionning
        - energies: dictionary to hold the results
        - dom_bounds: Boundaries of domains and PUs to calculate
    """
    j, domain = dom_bounds
    # Residues of the domain gathered little by little
    dom_residues = ""
    # Python multiprocessing is shit and cannot handle dict of dict.
    # So I have to use a tmp list for this little buddy...
    tmp_list = []
    for start_pu, end_pu in domain:
        pu_residues = ""
        dom_residues += ",".join([f"{str(x) + pdb_chain}" for x in range(start_pu, end_pu+1)]) + ","
        pu_residues += ",".join([f"{str(x) + pdb_chain}" for x in range(start_pu, end_pu+1)])
        pu_energy, pu_z_score = get_energy_and_z_score(f"{RESULTS_DIR}/{pdb_code_chain}", pu_residues)
        energies[(i, j, start_pu, end_pu)] = []
        tmp_list = energies[(i, j, start_pu, end_pu)]
        tmp_list.extend([pu_energy, pu_z_score])
        energies[(i, j, start_pu, end_pu)] = tmp_list
    dom_energy, dom_z_score = get_energy_and_z_score(f"{RESULTS_DIR}/{pdb_code_chain}", dom_residues)
    energies[(i, j)] = []
    tmp_list = energies[(i, j)]
    tmp_list.extend([dom_energy, dom_z_score])
    energies[(i, j)] = tmp_list


def write_peeling_results():
    """
    Parse Protein Peeling 3 results and calculate pseudo energy and AUL for
    all Protein Units.

    Args:
        - energies (dict): 
    """
    # Get old resnums
    peeling_num = os.path.join(RESULTS_DIR, "PDBs_Clean", pdb_code_chain, f"{pdb_code_chain}.num")
    with open(peeling_num, "r") as f:
        ori_resnums = [int(resnum) for resnum in f.readline().split()]
    # Parse Peeling.log
    peeling_results = {}
    peeling_log = os.path.join(RESULTS_DIR, "PDBs_Clean", pdb_code_chain, "Peeling", "Peeling.log")
    with open(peeling_log, "r") as f:
        # Index of the alternative partitionings
        nb_lvl = 1
        # Dict containing all the informations given by Peeling output
        for line in f:
            if not line.startswith("#") and line != "\n":
                line = line.split()
                peeling_results[nb_lvl] = {}
                peeling_results[nb_lvl]["i/e"] = float(line[0])
                peeling_results[nb_lvl]["i/i+e"] = float(line[1])
                peeling_results[nb_lvl]["R2"] = float(line[2])
                peeling_results[nb_lvl]["CI"] = float(line[3])
                peeling_results[nb_lvl]["N"] = int(line[4])
                # Retrieve only the PUs
                # store them as a list of tuples
                peeling_results[nb_lvl]["PUs"] = [(ori_resnums[int(line[5+i]) - 1], ori_resnums[int(line[5+i+1]) - 1]) for i in range(0, len(line[5:])-1, 2)]
                # Sorte by inceasing boundaries
                peeling_results[nb_lvl]["PUs"] = sorted(peeling_results[nb_lvl]["PUs"], key=lambda x: x[0]) 
                nb_lvl += 1
    
    # WRITE THE RESULTS
    logging.info("Calculate pseudo-energies of PUs")
    peeling = os.path.join(RESULTS_DIR, "PEELING_summary.txt")
    with open(peeling, "w") as f:
        for lvl, data in peeling_results.items():
            f.write(f"""Peeling level {lvl}\n    Number of Protein Units: {data["N"]}\n    Compaction Index: {round(data["CI"], 2)}\n""")
            for start_pu, end_pu in data["PUs"]:
                pu_residues = ",".join([f"{str(x) + pdb_chain}" for x in range(start_pu, end_pu+1)])
                pu_energy, pu_z_score = get_energy_and_z_score(f"{RESULTS_DIR}/{pdb_code_chain}", pu_residues)
                f.write(f"""    {str(start_pu)+"-"+str(end_pu):>7}: AUL={int((1-(1/(pu_z_score)**2))*100) if abs(pu_z_score) >= 1 else 0:3}% Z-score={round(pu_z_score, 1)}\n""")
    logging.info("Write the Peeling results")

if __name__ == '__main__':

    #################
    # Parse arguments
    #################

    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-u", "--uniprot-code", help="AlphaFold Uniprot Accession Id", type=str)
    group.add_argument("-p", "--pdb-code", help="PDB code", type=str)
    group.add_argument("-f", "--pdb-file", help="PDB file path", type=str)
    parser.add_argument("-c", "--pdb-chain", help="PDB chain. Default is A.", type=str, required=False, default="A")
    parser.add_argument("-o", "--output", help="Output directory", type=str, required=True)
    args = parser.parse_args()

    uniprot_code = args.uniprot_code
    pdb_code = args.pdb_code
    pdb_file = args.pdb_file
    pdb_chain = args.pdb_chain
    output_dir = args.output
    nb_cpu = multiprocessing.cpu_count()

    if pdb_file:
        if os.path.exists(pdb_file):
            pdb_code_chain = os.path.basename(os.path.splitext(pdb_file)[0]) + "_" + pdb_chain
        else:
            sys.exit("Unable to open file: " + pdb_file)
    elif uniprot_code:
        pdb_code_chain = uniprot_code + "_" + pdb_chain
    else:
        pdb_code_chain = pdb_code + "_" + pdb_chain

    confProDy(verbosity="none")

    # Configure the logger to redirect to job log file
    logging.basicConfig(level=logging.INFO,
                        handlers=[
                            logging.StreamHandler()
                        ],
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt="%Y/%m/%d %H:%M:%S")

    # Be conservative, if results directory already exists, create another one with suffix
    RESULTS_DIR = os.path.join(output_dir, pdb_code_chain)
    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR)
    else:
        name_rep = time.strftime("_%d_%m_%Y_") + ''.join(random.choice(list(map(str, list(range(0, 10))))) for i in range(5))
        RESULTS_DIR += name_rep
        os.makedirs(RESULTS_DIR)
        # Configure file log here
        fh = logging.FileHandler(os.path.join(RESULTS_DIR, "sword2.log"))
        fh.setLevel(logging.DEBUG)
        formatter = logging.Formatter("%(asctime)s %(levelname)s %(filename)s %(funcName)s() %(lineno)s - %(message)s")
        fh.setFormatter(formatter)
        logging.getLogger('').addHandler(fh)
        logging.warning(f"Results dir '{os.path.join(output_dir, pdb_code_chain)}' already exists. We created '{RESULTS_DIR}' instead.")


    ######################################
    # Define paths, variables, and logging
    ######################################

    BASE_DIR = os.path.abspath(os.path.dirname(__file__))
    BIN_DIR = os.path.join(BASE_DIR, "bin")
    SWORD_DIR = os.path.join(BIN_DIR, "SWORD/bin/SWORD")
    SWORD = os.path.join(SWORD_DIR, "SWORD")
    DSSP2PDB = os.path.join(BIN_DIR, "dssp2pdb.pl")
    DISPLAY_SWORD2 = os.path.join(BIN_DIR, "display_SWORD2_output.pl")

    # CHECK ENTRIES
    prot = check_parsing_pdb(uniprot_code, pdb_code, pdb_chain, pdb_file)

    # ESTIMATE THE RUNTIME
    est_time_in_minutes = int(predict_time(prot))
    if est_time_in_minutes < 1:
        est_time_in_minutes = 1

    logging.info(f"Estimated runtime: {est_time_in_minutes} minutes")

    ######################################################
    # Write the specific chain as a new PDB file for SWORD
    ######################################################

    pdb_chain_file = os.path.join(RESULTS_DIR, pdb_code_chain + ".pdb")
    logging.info("Write a clean version of the PDB: remove non standard residues")
    # Remove residues which have insertion codes
    res_to_remove = " "
    hv = prot.getHierView()
    for residue in hv.iterResidues():
        if residue.getIcode():
            res_to_remove += f"{residue.getResnum()} "
    # Keep only what constitutes the protein, remove the non standard residues
    if res_to_remove != " ":
        prot = prot.select('protein and not nonstdaa and not hetatm and not resnum ' + res_to_remove)
    else:
        prot = prot.select('protein and not nonstdaa and not hetatm')
    # In case there are negative residue numbers, we renumber
    # for convenience
    hv = prot.getHierView()
    SEQ_RESNUMS = [residue.getResnum() for residue in hv.iterResidues()]
    NEW_SEQ_RESNUMS = [i for i, _ in enumerate(SEQ_RESNUMS, start=1)]
    # Renumber the residues from 0 -> len(prot)
    [residue.setResnum(NEW_SEQ_RESNUMS[i]) for i, residue in enumerate(hv.iterResidues())]
    # Write the cleaned protein
    file_written = writePDB(pdb_chain_file, prot)
    if file_written is None:
        logging.warning(f"PDB file {pdb_chain_file} could not be written using function writePDB of ProDy")
        sys.exit(1)
    # Remove the ".pdb" extension for SWORD
    os.rename(file_written, os.path.splitext(file_written)[0])

    # Now that the protein is clean, get some infos
    hv = prot.getHierView()
    ch = hv.getChain(pdb_chain)
    # Get the aa sequence
    seq = ch.getSequence()
    PROT_LEN = len(seq)
    # Map the original resnums of the PDB to new numerotation from 0 -> N
    # this is used for the selection/coloration of domains and PUs in the sequence
    # section of the 3D viewer
    SEQ_RESNUMS_NORM = {num: i for i, num in enumerate(SEQ_RESNUMS)}
    SEQ_RESNUMS_NORM_INV = {i: num for i, num in enumerate(SEQ_RESNUMS)}

    ##############
    # Launch SWORD
    ##############

    logging.info("Launch SWORD")
    cmd_args = f"{DISPLAY_SWORD2} '{SWORD} -i {pdb_code_chain} --dir {RESULTS_DIR} -max 9'"
    cmd_args = shlex.split(cmd_args)
    output = subprocess.run(cmd_args, capture_output=True, check=True)
    output = output.stdout.decode("utf-8")
    output = output.split("\n")
    with open(f"{RESULTS_DIR}/sword.txt", "w") as f:
        for line in output:
            f.write(line + "\n")

    ####################
    # Parse SWORD output
    ####################

    sword_results = parse_sword(output)

    #####################################################
    # Calculate the energy and Z-score of PUs and Domains
    #####################################################

    logging.info("Calculate pseudo-energies of Domains")
    manager = multiprocessing.Manager()
    energies = manager.dict()
    for i, part in sword_results["DOMAINS"].items():
        with multiprocessing.Pool(processes=nb_cpu) as p:
            FUNC = partial(multiprocess_get_energy, i, energies)
            p.imap_unordered(FUNC, list(part["BOUNDARIES"].items()))
            p.close()
            p.join()

    ############################
    # Write  SWORD partitionings
    ############################

    write_partitionings(sword_results, energies)

    #######################
    # Write Peeling results
    #######################

    write_peeling_results()

    ############
    # Set colors
    ############

    pus_colors, dom_colors = define_colors(sword_results)

    ######################################################
    # Write Javascript formated histogram of SWORD Domains
    ######################################################

    write_domains_histogram(sword_results, dom_colors)

    #######################
    # Launch and parse DSSP
    #######################

    os.rename(os.path.splitext(pdb_chain_file)[0], pdb_chain_file)
    logging.info(f"Perform DSSP on {os.path.basename(pdb_chain_file)}")
    execDSSP(pdb_chain_file, outputname=f"{pdb_code_chain}", outputdir=RESULTS_DIR)
    prot = parseDSSP(os.path.join(RESULTS_DIR, pdb_code_chain + ".dssp"), parsePDB(pdb_chain_file, model=1))
    logging.info("dssp2pdb: translate DSSP output into a PDB HEADER") 
    # https://github.com/ivanamihalek/perlscr/blob/master/translation/dssp2pdb.pl
    dssp_output = os.path.join(RESULTS_DIR, f"{pdb_code_chain}.dssp")
    cmd_args = f"{DSSP2PDB} {dssp_output} {pdb_chain_file} > {pdb_chain_file}_dssp.pdb"
    cmd_args = os.system(cmd_args)

    #############
    # Contact map
    #############

    font = {'weight' : 'normal',
            'size'   : 11,}
    plt.rc('font', **font)
    plt.rcParams['axes.linewidth'] = 0.5
    plt.rcParams['xtick.major.size'] = 1.5
    plt.rcParams['ytick.major.size'] = 1.5
    plt.rcParams['figure.max_open_warning'] = 0

    os.makedirs(os.path.join(RESULTS_DIR, "Contact_Probability_Matrix"), exist_ok=True)
    # Load the contact map calculated by Peeling
    mat = np.loadtxt(os.path.join(RESULTS_DIR, "PDBs_Clean", pdb_code_chain, "file_proba_contact.mat"))
    for i, part in sword_results["DOMAINS"].items():
        # Draw the PUs of the alternative partition: 1 image per alternative
        adapt_height_image = 6.5 + (sum([len(domain) for j, domain in part["BOUNDARIES"].items()]) % 8) * 0.01
        fig1 = plt.figure(figsize=(5, adapt_height_image), dpi=300)
        ax1 = fig1.add_subplot()
        ax1.set_xlabel('Residues')
        ax1.set_ylabel('Residues')
        # Contact map
        ax1.imshow(mat, cmap='RdPu')
        # Invert Y axis to have both zeros at bottom left
        plt.gca().invert_yaxis()
        # Reduce padding on top of plot
        box1 = ax1.get_position()
        ax1.set_position([box1.x0, box1.y0 + box1.height * 0.2,
                        box1.width, box1.height * 0.9])
        if i == 0:
            ax1.set_title(f"Contact Probability Map of the\noptimal partition (all Protein Units)")
        else:
            ax1.set_title(f"Contact Probability Map of the alternative\npartition n°{i} (all Protein Units)")
        for j, domain in part["BOUNDARIES"].items():
            # Draw the PUs of domains: 1 image per domain
            adapt_height_image = 6.5 + (len(domain) % 8) * 0.02
            fig2 = plt.figure(figsize=(5, adapt_height_image), dpi=300)
            ax2 = fig2.add_subplot()
            ax2.set_xlabel('Residues')
            ax2.set_ylabel('Residues')
            # Contact map
            ax2.imshow(mat, cmap='RdPu')
            # Invert Y axis to have both zeros at bottom left
            plt.gca().invert_yaxis()
            box2 = ax1.get_position()
            ax2.set_position([box2.x0, box2.y0 + box2.height * 0.2,
                        box2.width, box2.height * 0.9])
            if i == 0:
                ax2.set_title(f"Contact Probability Map of the domain {j+1}\nof the optimal partition")
            else:
                ax2.set_title(f"Contact Probability Map of the domain {j+1}\nof the alternative partition n°{i}")
            for start_pu, end_pu in domain:
                # Save fig for individual PUs
                fig3 = plt.figure(figsize=(5, 6.5), dpi=300)
                ax3 = fig3.add_subplot()
                ax3.set_xlabel('Residues')
                ax3.set_ylabel('Residues')
                # Contact map
                ax3.imshow(mat, cmap='RdPu')
                plt.gca().invert_yaxis()
                ax3.set_title(f"Contact Probability Map of PU {start_pu}-{end_pu} of the domain {j+1}\nof the alternative partition n°{i}")
                # Draw PUs on contact map
                l = end_pu-start_pu
                rect = patches.Rectangle((start_pu-1, start_pu-1), l, l, linewidth=1.5, edgecolor="#%02x%02x%02x"%pus_colors[(start_pu, end_pu)], facecolor='none')
                rect.set_label(f"{start_pu}-{end_pu}")
                # need to copy the object to use it multiple times
                rect_copy = copy(rect)
                rect_copy2 = copy(rect)
                ax1.add_patch(rect)
                ax2.add_patch(rect_copy)
                ax3.add_patch(rect_copy2)
                ax3.legend(title="Protein Unit", loc='upper center', bbox_to_anchor=(0.5, -0.15),
                        fancybox=False, shadow=False, ncol=3)
                fig3.savefig(os.path.join(RESULTS_DIR, "Contact_Probability_Matrix", f"contact_probability_matrix_alternative_{i}_domain_{j}_pu_{start_pu}_{end_pu}.png"))
                plt.cla()
            # Put a legend below current axis
            ax2.legend(title="Protein Units", loc='upper center', bbox_to_anchor=(0.5, -0.15),
                    fancybox=True, shadow=False, ncol=3)
            fig2.savefig(os.path.join(RESULTS_DIR, "Contact_Probability_Matrix", f"contact_probability_matrix_alternative_{i}_domain_{j}.png"))
            plt.cla()
        # Put a legend below current axis
        ax1.legend(title="Protein Units", loc='upper center', bbox_to_anchor=(0.5, -0.15),
                fancybox=True, shadow=False, ncol=3)
        fig1.savefig(os.path.join(RESULTS_DIR, "Contact_Probability_Matrix", f"contact_probability_matrix_alternative_{i}.png"))
        plt.cla()


    #########################
    # Junctions consistencies
    #########################

    logging.info("Calculate junctions consistencies")

    cmd_args = f"{BIN_DIR}/stat_pu_domains_from_SWORD.pl {RESULTS_DIR}/sword.txt"
    cmd_args = shlex.split(cmd_args)
    output = subprocess.run(cmd_args, capture_output=True, check=True)
    output = output.stdout.decode("utf-8")
    output = output.split("\n")
    junctions = {}
    with open(f"{RESULTS_DIR}/junctions_consistencies.txt", "w") as f:
        for line in output:
            jctn_found = re.search(r"^(\d+)\s+(\d+)\s+(\d\.\d+)\s+(\d+\.\d+)$", line)
            if jctn_found:
                jct = int(jctn_found.group(1))
                cnt = int(jctn_found.group(2))
                raw = float(jctn_found.group(3))
                weight = float(jctn_found.group(4))
                junctions[jct] = { "cnt": cnt, "raw": raw, "weight": weight }
            f.write(line + "\n")
        # Remove the first and the last junctions
        # because they simply are the beginning and the end of the sequence
        # so not relevant to show
        del junctions[min(list(junctions.keys()))]
        del junctions[max(list(junctions.keys()))]


    # Mapping of authors PDB residues numbers with the new one presented on the server 
    with open(f"{RESULTS_DIR}/mapping_auth_resnums.txt", "w") as f1:
        f1.write("#Mapping of authors PDB residues numbers with the new one presented on the server\nORIGINAL RENUM\n")
        for i, j in enumerate(SEQ_RESNUMS, start=1):
            f1.write(f"{j} {i}\n")

    # Clean the directory
    logging.info("Clean and prepare results")
    shutil.rmtree(os.path.join(RESULTS_DIR, "PDBs_Stand"), ignore_errors=False, onerror=None)
    shutil.move(os.path.join(RESULTS_DIR, "PDBs_Clean"), os.path.join(RESULTS_DIR, "SWORD"))
    os.makedirs(os.path.join(RESULTS_DIR, "Junctions"), exist_ok=True)
    shutil.move(os.path.join(RESULTS_DIR, "junctions_consistencies.txt"), os.path.join(RESULTS_DIR, "Junctions"))
    os.system(f"mv {RESULTS_DIR}/SWORD/*/Peeling {RESULTS_DIR}/Protein_Units")

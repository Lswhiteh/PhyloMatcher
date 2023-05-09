## GBIF MATCHER

import argparse
import multiprocessing as mp
from tqdm import tqdm
from pygbif import species
import os


def gbif_read_csv(csvfile):
    target_list = []
    with open(csvfile, "r") as ifile:
        for line in ifile.readlines():
            if "Tree_Sp_Name" in line:
                continue
            else:
                target_list.append(line.strip().split(",")[0])

    return target_list


def get_sp_id(sp):
    sp_dict = species.name_backbone(sp)
    if "speciesKey" in sp_dict:
        return sp_dict["speciesKey"], sp_dict["species"]
    else:
        return None, None


def get_synonyms(sp_key):
    nu = species.name_usage(sp_key, data="synonyms")["results"]
    if len(nu) > 0:
        return [i["canonicalName"] for i in nu if "canonicalName" in i]
    else:
        return []


def worker(sp):
    key, curr_name = get_sp_id(sp)
    if key:
        synonyms = get_synonyms(key)
        synonyms.insert(0, sp)
        synonyms.append(curr_name)
    else:
        synonyms = [sp]

    return synonyms


def gbif_main(input_csv, outfile, threads):
    sp_list = gbif_read_csv(input_csv)
    cleaned_sp_list = [i.replace("_", " ") for i in sp_list]

    with mp.Pool(threads) as p:
        synonyms = list(
            tqdm(
                p.imap(worker, cleaned_sp_list, chunksize=4),
                desc="[INFO] Fetching GBIF information",
                total=len(cleaned_sp_list),
            )
        )

    max_len = max([len(i) for i in synonyms])
    eq_headers = (
        ["Tree_Sp_Name"] + [f"Eq_{i}" for i in range(max_len - 1)] + ["Curr_Name"]
    )

    os.makedirs(os.path.dirname(outfile), exist_ok=True)

    with open(outfile, "w") as ofile:
        ofile.write("\t".join(eq_headers) + "\n")
        for names in synonyms:
            ofile.write("\t".join([i.replace(" ", "_") for i in names]) + "\n")


## ENTREZ MATCHER

import argparse
import os
import Bio

from tqdm import tqdm
from Bio import Entrez


def get_ua():
    ap = argparse.ArgumentParser()
    ap.add_argument("-e", "--email", dest="email", required=True)
    ap.add_argument(
        "-c",
        "--csv",
        dest="input_csv",
        required=True,
        help="CSV where first column is a list of target species names to look up. Will check other columns too. Ignores final column, assumes is notes.",
    )
    return ap.parse_args()


def run_esearch(sp):
    try:
        e_handle = Entrez.esearch(db="taxonomy", term=sp, rettype="gb")
        record = Entrez.read(e_handle)["IdList"]
        e_handle.close()
        return record[0]
    except:
        raise Exception


def get_tax_id(species):
    """Searches Entrez Taxonomy and grabs the official ID"""
    try:
        sp_id = run_esearch(species)
        return sp_id
    except:
        # Check for misspellings
        if species[-2:] == "um":
            species = "us".join(species.rsplit("um", 1))
        elif species[-2:] == "us":
            species = "um".join(species.rsplit("us", 1))
        elif species[-2:] == "ea":
            species = "eus".join(species.rsplit("ea", 1))
        elif species[-2:] == "eus":
            species = "ea".join(species.rsplit("eus", 1))
        try:
            sp_id = run_esearch(species)
            return sp_id
        except:
            return None


def get_tax_info(tax_ids):
    fetch = Entrez.efetch(id=tax_ids, db="taxonomy")
    taxinfo = Entrez.read(fetch)
    return taxinfo


def get_other_names(tax_dict):
    synonyms = []
    if "OtherNames" in tax_dict:
        if "Synonym" in tax_dict["OtherNames"]:
            synonyms.extend(tax_dict["OtherNames"]["Synonym"])

        if "Name" in tax_dict["OtherNames"]:
            for name in tax_dict["OtherNames"]["Name"]:
                if name["ClassCDE"] == "authority":
                    synonyms.append(" ".join(name["DispName"].split()[:2]))

    return list(set(synonyms))


def ncbi_read_csv(csvfile):
    sp_list = []
    target_list = []
    idx_list = []
    with open(csvfile, "r") as ifile:
        for idx, line in enumerate(ifile.readlines(), start=-1):
            if "Tree_Sp_Name" in line:
                continue
            else:
                _line = [i for i in line.strip().split(",")[:-1] if i]
                sp_list.extend(_line)
                target_list.append(line.strip().split(",")[0])
                idx_list.extend([idx] * len(_line))

    return sp_list, target_list, idx_list


def entrez_main(input_csv, outfile, user_email):
    Entrez.email = user_email
    os.makedirs(os.path.dirname(outfile), exist_ok=True)

    sp_list, target_list, idx_list = ncbi_read_csv(input_csv)
    cleaned_sp_list = [i.replace("_", " ") for i in sp_list]

    name_res = []
    pass_list = []
    fail_list = []
    tax_ids = []
    for idx, sp_str in tqdm(
        enumerate(cleaned_sp_list), "[INFO] Fetching IDs", total=len(cleaned_sp_list)
    ):
        try:
            tax_id = get_tax_id(sp_str)
            if tax_id:
                tax_ids.append(get_tax_id(sp_str))
                pass_list.append(sp_str)
            else:
                fail_list.append(target_list[idx_list[idx]])

        except:
            fail_list.append(target_list[idx_list[idx]])

    # Filter out redundant hits

    unique_tax_ids = set(tax_ids)
    pass_list = [pass_list[i] for i in [tax_ids.index(j) for j in unique_tax_ids]]

    tax_info = get_tax_info(unique_tax_ids)

    for name, ti in tqdm(
        zip(pass_list, tax_info), "[INFO] Parsing XML data", total=len(pass_list)
    ):
        try:
            namelist = get_other_names(ti)
        except:
            fail_list.append(name)
            continue
        if name in namelist:
            namelist.remove(name)

        namelist = [i.replace(" ", "_") for i in namelist]
        namelist.append(name.replace(" ", "_"))
        name_res.append(namelist)

    # Remove redundant lists to allow for checking all csv entries
    name_res = set([tuple(i) for i in name_res])
    name_res = [list(i) for i in name_res]

    for n in name_res:
        for s in target_list:
            if s in n:
                n.remove(s)
                n.insert(0, s)

    max_len = max([len(i) for i in name_res])
    eq_headers = ["Tree_Sp_Name"] + [f"Eq_{i}" for i in range(max_len - 1)]

    with open(outfile, "w") as ofile:
        ofile.write("\t".join(eq_headers) + "\n")
        for names in name_res:
            ofile.write("\t".join(names) + "\n")

    with open(f"{outfile.split('.')[0]}_fails.csv", "w") as failfile:
        for i in set(fail_list):
            failfile.write(i.replace(" ", "_") + "\n")

    print("[INFO] Done")


## TRAITMATCHER

import pandas as pd
import csv
from tqdm import tqdm
import os


def trait_main(traitfile, speciesfile, outfile):
    os.makedirs(os.path.dirname(outfile), exist_ok=True)

    with open(traitfile, "r") as tfile:
        dialect = csv.Sniffer().sniff(tfile.readline())

    trait_df = pd.read_csv(traitfile, delimiter=dialect.delimiter)
    species_df = pd.read_csv(speciesfile)

    specieslist = set(list(species_df.stack().values))
    for sp in tqdm(specieslist, desc="Matching taxa to traits"):
        trait_df.loc[trait_df.iloc[:, 0] == sp, "tree_name"] = sp

    trait_df.to_csv(outfile, index=False)


## GUI

import PySimpleGUI as sg
import os
import shutil

sg.theme("DarkAmber")  # set GUI color

# Define all the stuff inside GUI window
layout = [
    [sg.Text("Enter CSV file:")],
    [
        sg.Input(enable_events=True, key="-IN_CSV-"),
        sg.FileBrowse(),
    ],  # Let GUI detect inputs and add browse button
    [sg.Text("Select a destination folder:", visible=False, key="-OUTPUT_MSG-")],
    [
        sg.Input(enable_events=True, visible=False, key="-DEST-"),
        sg.FolderBrowse(visible=False, key="-DEST_BROWSE-"),
    ],
    [sg.Text("", key="options")],  # To reveal options later
    [sg.Button("Run Taxomatcher", visible=False, enable_events=True)],
    [sg.Text("", key="options_2", size=(1))],
    [sg.Button("Run Traitmatcher", visible=False, enable_events=True)],
    [sg.Button("Help"), sg.Button("License"), sg.Button("Exit")],
]

# Create the GUI window
window = sg.Window("Taxomatcher", layout, size=(550, 450))

config = True

# Event Loop to process GUI events
while True:
    event, values = window.read()
    if (
        event == sg.WIN_CLOSED or event == "Exit"
    ):  # if user closes window or clicks exit
        break

    if event == "License":
        license_text = """MIT License

            Copyright (c) 2023 Logan Whitehouse

            Permission is hereby granted, free of charge, to any person obtaining a copy
            of this software and associated documentation files (the "Software"), to deal
            in the Software without restriction, including without limitation the rights
            to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
            copies of the Software, and to permit persons to whom the Software is
            furnished to do so, subject to the following conditions:

            The above copyright notice and this permission notice shall be included in all
            copies or substantial portions of the Software.

            THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
            IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
            FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
            AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
            LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
            OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
            SOFTWARE."""

        license_window = sg.Window("License", [[sg.Text(license_text)]])
        while True:
            event, values = license_window.read()
            if event in (None, "OK"):
                break
        license_window.close()

    # If press help
    if event == "Help":
        help_layout = [
            [
                sg.TabGroup(
                    [
                        [
                            sg.Tab(
                                "Main",
                                [
                                    [
                                        sg.Multiline(
                                            "Make sure none of the files being handled by the program are currently open (ex: in excel).\n\nSpecies name input file should be single column CSV of target species to look up synonyms for (other columns will be ignored). Names can be either space or underscore ('_') seperated - 'Sphenodon_punctatus' is equivalent to 'Sphenodon punctatus'. See 'Example Sp_CSV' tab.\n\nSpecies trait input file should be CSV with species name in first column. Species names should be seperated by underscores and case-sensetive. See 'Example Trait_CSV' tab.\n\nFor more info, see: https://github.com/Lswhiteh/taxomatcher/blob/main/README.md",
                                            size=(70, 10),
                                            key="-MAIN_TEXT-",
                                        )
                                    ],
                                    [sg.Button("OK")],
                                ],
                            )
                        ],
                        [
                            sg.Tab(
                                "Example Sp_CSV",
                                [
                                    [
                                        sg.Table(
                                            values=[
                                                ["Homo sapiens"],
                                                ["Gecko gekko"],
                                                ["Gorilla gorilla"],
                                            ],
                                            headings=["Species Name"],
                                            auto_size_columns=True,
                                            num_rows=10,
                                            vertical_scroll_only=True,
                                        )
                                    ]
                                ],
                            )
                        ],
                        [
                            sg.Tab(
                                "Example Trait_CSV",
                                [
                                    [
                                        sg.Table(
                                            values=[
                                                ["Homo_sapiens", "Omnivore", "Diurnal"],
                                                [
                                                    "Gecko_gekko",
                                                    "Carnivore",
                                                    "Nocturnal",
                                                ],
                                                [
                                                    "Gorilla_gorilla",
                                                    "Omnivore",
                                                    "Diurnal",
                                                ],
                                            ],
                                            headings=[
                                                "Species Name",
                                                "Diet",
                                                "Activity",
                                            ],
                                            auto_size_columns=True,
                                            num_rows=10,
                                            vertical_scroll_only=True,
                                        )
                                    ]
                                ],
                            )
                        ],
                    ]
                )
            ]
        ]
        help_window = sg.Window("Help", help_layout)
        while True:
            event, values = help_window.read()
            if event in (None, "OK"):
                break
        help_window.close()

    if "-IN_CSV-" in values:
        if values["-IN_CSV-"]:
            # Prompt the user to select a file destination
            window["-OUTPUT_MSG-"].update(visible=True)
            window["-DEST-"].update(visible=True)
            window["-DEST_BROWSE-"].update(visible=True)

    # If they've selected input file & destination, show options + run button
    if "-IN_CSV-" in values and "-DEST-" in values and config == True:
        search_layout = [
            [
                sg.Radio(
                    "Search GBIF (recommended)", "option", key="gbif", default=True
                ),
                sg.Radio(
                    "Search NCBI (many misses, not recommended)", "option", key="ncbi"
                ),
            ],
            [
                sg.Text("Enter your professional email:"),
                sg.Input(enable_events=True, size=(25, 1), key="-EMAIL-"),
            ],
            [
                sg.Text("Enter a thread count (optional):"),
                sg.Input(enable_events=True, size=(5, 1), key="-THREADS-"),
            ],
        ]
        window.extend_layout(window["options"], search_layout)

        # Add run button
        window["Run Taxomatcher"].update(visible=True)

        window.finalize()
        config = False

    # If they click run
    if event == "Run Taxomatcher":
        # Check if a file has been selected:
        if not values["-IN_CSV-"] or not values["-DEST-"]:
            sg.popup("Please select a file to process and a place to store the output.")
            continue

        elif not values["-IN_CSV-"].endswith(".csv"):
            sg.popup(
                'Taxomatcher only supports analysis of CSV files. Please input a CSV, and see "Help" for more details.'
            )
            continue

        else:
            chosen_csv = values["-IN_CSV-"]
            dest_dr = values["-DEST-"]
            run_name = chosen_csv.split("/")[-1].split(".")[0]

            # Make final output file destination from given dir
            if values["gbif"]:
                taxo_file = f"{values['-DEST-']}/{run_name}_GBIF_output.tsv"

            elif values["ncbi"]:
                taxo_file = f"{values['-DEST-']}/{run_name}_NCBI_output.tsv"

            # Use their chosen thread value, else default to 4
            if values["-THREADS-"] != "":
                chosen_threads = values["-THREADS-"]
            else:
                chosen_threads = 4

            lock_configs = {
                window["gbif"].update(disabled=True),
                window["ncbi"].update(disabled=True),
                window["-EMAIL-"].update(disabled=True),
                window["-THREADS-"].update(disabled=True),
                window["Run Taxomatcher"].update(disabled=True),
            }

            # Run Taxomatcher on their chosen database
            if values["gbif"]:
                lock_configs
                gbif_main(
                    input_csv=chosen_csv, outfile=taxo_file, threads=chosen_threads
                )

            elif values["ncbi"]:
                if values["-EMAIL-"] == "":
                    sg.popup("NCBI requires entry of a user email to query.")
                    continue
                else:
                    lock_configs
                    entrez_main(
                        input_csv=chosen_csv,
                        outfile=taxo_file,
                        user_email=values["-EMAIL-"],
                    )

            sg.popup(
                f"Taxomatcher has completed successfully. Results have been saved to {dest_dr}"
            )

            trait_layout = [
                [
                    sg.Text(
                        "Optional: to run trait matcher, enter trait file below (see Help for more info)"
                    )
                ],
                [sg.Input(enable_events=True, key="-TRAIT_FILE-"), sg.FileBrowse()],
                [sg.Button("Run Traitmatcher", visible=False, enable_events=True)],
            ]
            window.extend_layout(window["options_2"], trait_layout)
            window.finalize()

    if "-TRAIT_FILE-" in values:
        if values["-TRAIT_FILE-"]:
            window["Run Traitmatcher"].update(visible=True)

    if event == "Run Traitmatcher":
        # Save trait file in same place as saved taxo output
        trait_file = f"{values['-DEST-']}/{run_name}_traitmatch_output.csv"

        trait_main(
            traitfile=values["-TRAIT_FILE-"], speciesfile=taxo_file, outfile=trait_file
        )

        sg.popup(
            f"Traitmatcher has completed successfully. Results have been saved to {dest_dr} (same as Taxomatcher output)"
        )
        break

window.close()

import PySimpleGUI as sg
import os

from gbif_matcher import main as gbif_main
from entrez_matcher import main as entrez_main
from trait_matcher import main as trait_main

sg.theme("DarkAmber")  

# Define PhyloMatcher GUI
taxo_layout = [
    [sg.Text("Enter CSV file:")],
    [
        sg.Input(enable_events=True, key="-IN_CSV-"),
        sg.FileBrowse(),
    ],
    [sg.Text("Select a destination folder:", visible=False, key="-OUTPUT_MSG-")],
    [
        sg.Input(enable_events=True, visible=False, key="-DEST-"),
        sg.FolderBrowse(visible=False, key="-DEST_BROWSE-"),
    ],
    [sg.Text("", key="options")],  # To reveal taxomatch options later
    [sg.Button("Run PhyloMatcher", visible=False, enable_events=True)],
    [sg.Text("", key="options_2", size=(1))],
    [sg.Button("Run Traitmatcher", visible=False, enable_events=True)],
    [sg.Button("Help"), sg.Button("License"), sg.Button("Exit")],
]

# Define standalone traitmatcher GUI
trait_layout = [
    [sg.Text("Enter a trait CSV file:")],
    [
        sg.Input(enable_events=True, key="-IN_CSV-"),
        sg.FileBrowse(),
    ],

    [sg.Text("Enter a species TSV file:")],
    [
        sg.Input(enable_events=True, key="-SPEC_FILE-"),
        sg.FileBrowse(key="in_spec_file"),
    ],

    [sg.Text("Select a destination folder:")],
    [
        sg.Input(enable_events=True, visible=True, key="-DEST-"),
        sg.FolderBrowse(visible=True, key="-DEST_BROWSE-"),
    ],

    [sg.Button("Run Traitmatcher (standalone)", visible=False, enable_events=True)],
    [sg.Text("", key="trait_options")],  # To reveal standalone traitmatcher options later
]

# GUI to choose between PhyloMatcher and standalone traitmatcher
radio_layout = [
    [sg.Text("Choose a tool:")],
    [sg.Radio("PhyloMatcher", "RADIO", default=True, key="-TAXO-")],
    [sg.Radio("Standalone Traitmatcher", "RADIO", key="-TRAIT-")],
    [sg.Button("OK")]
]

# Create the initial GUI window
window = sg.Window("PhyloMatcher", radio_layout, size=(250, 150))

current_window = "choose"
config = True
trait_warned = False

# Event Loop to process GUI events
while True:
    event, values = window.read()
    if (
        event == sg.WIN_CLOSED or event == "Exit"
    ):  # if user closes window or clicks exit
        break

    # If the user clicks OK on the radio button window, close the radio window and show the chosen tool GUI
    if event == "OK":
        window.close()

        # If PhyloMatcher is chosen, show PhyloMatcher GUI
        if values["-TAXO-"]:
            window.close()
            window = sg.Window("PhyloMatcher", taxo_layout, size=(550, 450))
            current_window = "taxo"

        # If standalone traitmatcher is chosen, show standalone traitmatcher GUI
        elif values["-TRAIT-"]:
            window.close()
            window = sg.Window("Traitmatcher (standalone)", trait_layout, size=(550, 450))
            current_window = "trait"

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
                                            "Make sure none of the files being handled by the program are currently open (ex: in excel).\n\nSpecies name input file should be single column CSV of target species to look up synonyms for (other columns will be ignored). Names can be either space or underscore ('_') seperated - 'Sphenodon_punctatus' is equivalent to 'Sphenodon punctatus'. See 'Example Sp_CSV' tab.\n\nSpecies trait input file should be CSV with species name in first column. Species names should be seperated by underscores and case-sensetive. See 'Example Trait_CSV' tab.\n\nFor more info, see: https://github.com/Lswhiteh/PhyloMatcher/blob/main/README.md",
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
                        [
                            sg.Tab(
                                "Example TSV for standalone Traitmatcher",
                                [
                                    [
                                        sg.Table(
                                            values=[
                                                ["Homo_sapiens", "Homo_drennani", "Homo_columbicus"],
                                                [
                                                    "Gecko_gekko",
                                                ],
                                                [
                                                    "Gorilla_gorilla",
                                                    "Troglodytes_gorilla",
                                                ],
                                            ],
                                            headings=[
                                                "'Canonical' name",
                                                "Synonym 1",
                                                "Synonym 2",
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

    if current_window == "trait":
        if values.get("-SPEC_FILE-") is not None and trait_warned == False:
            sg.popup("Note: the species file must be a TSV formatted identical to PhyloMatcher output ('canonical'/tree names in first column, synonyms in subsequent rows; case-sensitive and separated by underscores) - see Help for example.")
            trait_warned = True
        
        # If all reqs filled, show traitmatcher
        required_inputs = ["-SPEC_FILE-", "-IN_CSV-", "-DEST-"]
        if all(key in values and values[key] for key in required_inputs):
            window["Run Traitmatcher (standalone)"].update(visible=True)
        
        if event == "Run Traitmatcher (standalone)":
            if not values["-SPEC_FILE-"].endswith(".tsv"):
                sg.popup(
                    'Traitmatcher only supports analysis of TSV files. Please input a TSV, and see "Help" for more details.'
                )
                continue
            else: 
                # Save trait file to dest dr
                run_name = values['-IN_CSV-'].split("/")[-1].split(".")[0]
                trait_file = f"{values['-DEST-']}/{run_name}_traitmatch_output.csv"

            trait_main(
                traitfile=values["-IN_CSV-"], speciesfile=values["-SPEC_FILE-"], outfile=trait_file
            )

            sg.popup(
                f"Traitmatcher has completed successfully. Results have been saved to {values['-DEST-']}."
            )
            break



    if current_window == "taxo":
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
            window["Run PhyloMatcher"].update(visible=True)

            window.finalize()
            config = False

        # If they click run
        if event == "Run PhyloMatcher":
            # Check if a file has been selected:
            if not values["-IN_CSV-"] or not values["-DEST-"]:
                sg.popup("Please select a file to process and a place to store the output.")
                continue

            elif not values["-IN_CSV-"].endswith(".csv"):
                sg.popup(
                    'PhyloMatcher only supports analysis of CSV files. Please input a CSV, and see "Help" for more details.'
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

                # Run PhyloMatcher on their chosen database
                if values["gbif"]:
                    lock_configs = {
                        window["gbif"].update(disabled=True),
                        window["ncbi"].update(disabled=True),
                        window["-EMAIL-"].update(disabled=True),
                        window["-THREADS-"].update(disabled=True),
                        window["Run PhyloMatcher"].update(disabled=True),
                    }
                    lock_configs  
                    gbif_main(
                        input_csv=chosen_csv, outfile=taxo_file, threads=chosen_threads
                    )

                elif values["ncbi"]:
                    if values["-EMAIL-"] == "":
                        sg.popup("NCBI requires entry of a user email to query.")
                        continue
                    else:
                        lock_configs = {
                            window["gbif"].update(disabled=True),
                            window["ncbi"].update(disabled=True),
                            window["-EMAIL-"].update(disabled=True),
                            window["-THREADS-"].update(disabled=True),
                            window["Run PhyloMatcher"].update(disabled=True),
                        }
                        lock_configs # defining twice is gross but it fixes annoying bug :(
                        entrez_main(
                            input_csv=chosen_csv,
                            outfile=taxo_file,
                            user_email=values["-EMAIL-"],
                        )

                sg.popup(
                    f"PhyloMatcher has completed successfully. Results have been saved to {dest_dr}"
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
                f"Traitmatcher has completed successfully. Results have been saved to {dest_dr} (same as PhyloMatcher output)"
            )
            break

window.close()
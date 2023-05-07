import PySimpleGUI as sg
import os
import shutil
from gbif_matcher import main as gbif_main
from entrez_matcher import main as entrez_main
from trait_matcher import main as trait_main

sg.theme('DarkAmber')   # set GUI color

# Define all the stuff inside GUI window
layout = [[sg.Text('Enter CSV file:')],
          
          [sg.Input(enable_events=True, key='-IN_CSV-'), sg.FileBrowse()], # Let GUI detect inputs and add browse button
    
          [sg.Text('Select a destination folder:', visible=False, key='-OUTPUT_MSG-')],
          [sg.Input(enable_events=True, visible=False, key='-DEST-'), sg.FolderBrowse(visible=False, key='-DEST_BROWSE-')],

          [sg.Text('', key='options')], # To reveal options later
          [sg.Button('Run Taxomatcher', visible=False, enable_events=True)],

          [sg.Text('', key='options_2', size=(1))],
          [sg.Button('Run Traitmatcher', visible=False, enable_events=True)],

          [sg.Button('Help'), sg.Button('Exit')],
          ]

# Create the GUI window
window = sg.Window('Taxomatcher', layout, size = (550,450))

config = True


# Event Loop to process GUI events
while True:
    event, values = window.read()
    if event == sg.WIN_CLOSED or event == 'Exit': # if user closes window or clicks exit
        break

    # If they press help
    if event == 'Help':
        help_layout = [
            [sg.TabGroup([
                [sg.Tab('Main', [
                    [sg.Multiline("Make sure none of the files being handled by the program are currently open (ex: in excel).\n\nSpecies name input file should be single column CSV of target species to look up synonyms for (other columns will be ignored). Names can be either space or underscore ('_') seperated - 'Sphenodon_punctatus' is equivalent to 'Sphenodon punctatus'. See 'Example Sp_CSV' tab.\n\nSpecies trait input file should be CSV with species name in first column. Species names should be seperated by underscores and case-sensetive. See 'Example Trait_CSV' tab.\n\nFor more info, see: https://github.com/Lswhiteh/taxomatcher/blob/main/README.md", size=(70, 10), key='-MAIN_TEXT-')],
                    [sg.Button('OK')]
                ])],
                [sg.Tab('Example Sp_CSV', [
                    [sg.Table(
                        values=[['Homo sapiens'], ['Gecko gekko'], ['Gorilla gorilla']],
                        headings=['Species Name'],
                        auto_size_columns=True,
                        num_rows=10,
                        vertical_scroll_only=True,
                    )]
                ])],
                [sg.Tab('Example Trait_CSV', [
                    [sg.Table(
                        values=[['Homo_sapiens', 'Omnivore', 'Diurnal'], 
                                ['Gecko_gekko', 'Carnivore', 'Nocturnal'], 
                                ['Gorilla_gorilla', 'Omnivore', 'Diurnal']],
                        headings=['Species Name', 'Diet', 'Activity'],
                        auto_size_columns=True,
                        num_rows=10,
                        vertical_scroll_only=True,
                    )]
                ])]
            ])]]
        help_window = sg.Window('Help', help_layout)
        while True:
            event, values = help_window.read()
            if event in (None, 'OK'):
                break
        help_window.close()


    if '-IN_CSV-' in values: 
        if values['-IN_CSV-']:
            # Prompt the user to select a file destination
            window['-OUTPUT_MSG-'].update(visible=True)
            window['-DEST-'].update(visible=True)
            window['-DEST_BROWSE-'].update(visible=True)

    # If they've selected input file & destination, show options + run button
    if '-IN_CSV-' in values and '-DEST-' in values and config == True:
            search_layout = [[sg.Radio('Search GBIF (recommended)', 'option', key='gbif', default=True),
                              sg.Radio('Search NCBI (many misses, not recommended)', 'option', key='ncbi')],
                             [sg.Text('Enter your professional email:'),
                              sg.Input(enable_events=True, size=(25, 1), key='-EMAIL-')],
                             [sg.Text('Enter a thread count (optional):'),
                              sg.Input(enable_events=True, size=(5, 1), key='-THREADS-')],
            ]
            window.extend_layout(window['options'], search_layout)

            # Add run button
            window['Run Taxomatcher'].update(visible=True)
            
            window.finalize()
            config = False

    # If they click run
    if event == 'Run Taxomatcher':
        # Check if a file has been selected:
        if not values['-IN_CSV-'] and not values['-DEST-']:
            sg.popup('Please select a file to process and a place to store the output.')
            continue
        
        elif not values['-IN_CSV-'].endswith('.csv'):
            sg.popup('Taxomatcher only supports analysis of CSV files. Please input a CSV, and see "Help" for more details.')
            continue
                
        else:
            chosen_csv = values['-IN_CSV-']
            dest_dr = values['-DEST-']
            run_name = chosen_csv.split("/")[-1].split(".")[0]

            # Make final output file destination from given dir
            if values['gbif']:
                taxo_file = f"{values['-DEST-']}/{run_name}_GBIF_output.tsv"
                
            elif values['ncbi']:
                taxo_file = f"{values['-DEST-']}/{run_name}_NCBI_output.tsv"

            # Use their chosen thread value, else default to 4
            if values['-THREADS-'] != '':
                chosen_threads = values['-THREADS-']
            else: 
                chosen_threads = 4


            lock_configs = {
                window['gbif'].update(disabled=True), 
                window['ncbi'].update(disabled=True), 
                window['-EMAIL-'].update(disabled=True), 
                window['-THREADS-'].update(disabled=True), 
                window['Run Taxomatcher'].update(disabled=True)}

            # Run Taxomatcher on their chosen database
            if values['gbif']:
                lock_configs
                gbif_main(input_csv = chosen_csv, outfile = taxo_file, threads = chosen_threads)

            elif values['ncbi']:
                if values['-EMAIL-'] == '':
                    sg.popup('NCBI requires entry of a user email to query.')
                    continue
                else: 
                    lock_configs
                    entrez_main(input_csv = chosen_csv, outfile = taxo_file, user_email=values['-EMAIL-'])

            sg.popup(f'Taxomatcher has completed successfully. Results have been saved to {dest_dr}')
            
            trait_layout = [ 
                             [sg.Text('Optional: to run trait matcher, enter trait file below (see Help for more info)')],
                             [sg.Input(enable_events=True, key='-TRAIT_FILE-'), sg.FileBrowse()],
                             [sg.Button('Run Traitmatcher', visible=False, enable_events=True)],
            ]
            window.extend_layout(window['options_2'], trait_layout)
            window.finalize()
        
    if '-TRAIT_FILE-' in values:
        if values['-TRAIT_FILE-']:
            window['Run Traitmatcher'].update(visible=True)
    
    if event == 'Run Traitmatcher':
        # Save trait file in same place as saved taxo output
        trait_file = f"{values['-DEST-']}/{run_name}_traitmatch_output.csv"

        trait_main(traitfile = values['-TRAIT_FILE-'], speciesfile = taxo_file, outfile = trait_file)

        sg.popup(f'Traitmatcher has completed successfully. Results have been saved to {dest_dr} (same as Taxomatcher output)')
        break 









window.close()
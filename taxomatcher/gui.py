import PySimpleGUI as sg
import os
import shutil
from gbif_matcher import main as gbif_main
from entrez_matcher import main as entrez_main

sg.theme('DarkAmber')   # set GUI color

# Define all the stuff inside GUI window
layout = [[sg.Text('Enter CSV file:')],
          
          [sg.Input(enable_events=True, key='-IN_CSV-'), sg.FileBrowse()], # Let GUI detect inputs and add browse button
          
          [sg.Text('', key='options')], # To reveal options later
          [sg.Button('Run Taxomatcher', visible=False, enable_events=True)],
    
          [sg.Text('', size=(1))], # Add a line space
          [sg.Text('Select a destination folder:', visible=False, key='-OUTPUT_MSG-')],
          [sg.Input(enable_events=True, visible=False, key='-DEST-'), sg.FolderBrowse(visible=False, key='-DEST_BROWSE-')],

          [sg.Text('', size=(1))],
          [sg.Button('Save', visible=False), sg.Button('Help'), sg.Button('Exit')],
          ]

# Create the GUI window
window = sg.Window('Taxomatcher', layout, size = (550,350))

config = True

# Event Loop to process GUI events
while True:
    event, values = window.read()
    if event == sg.WIN_CLOSED or event == 'Exit': # if user closes window or clicks exit
        break

    # If they press help
    if event == 'Help':
        sg.popup('CSV input file should be single column of target species to look up (other columns will be ignored). Names can be either space or underscore ("_") separated, i.e. "Sphenodon_punctatus" is equivalent to "Sphenodon punctatus". For more info, see https://github.com/Lswhiteh/taxo-matcher')

    # If they've selected input file, show config options
    if values['-IN_CSV-'] and config == True:
            search_layout = [[sg.Radio('Search GBIF (recommended)', 'option', key='gbif', default=True),
                              sg.Radio('Search NCBI (many misses, not recommended)', 'option', key='ncbi')],
                             [sg.Text('Enter your professional email:'),
                              sg.Input(enable_events=True, size=(25, 1), key='-EMAIL-')],
                             [sg.Text('Enter a thread count (optional):'),
                              sg.Input(enable_events=True, size=(5, 1), key='-THREADS-')],
            ]
            window.extend_layout(window['options'], search_layout)
            window['Run Taxomatcher'].update(visible=True)
            window.finalize()
            config = False

    # If they click run
    if event == 'Run Taxomatcher':
        # Check if a file has been selected:
        if not values['-IN_CSV-']:
            sg.popup('Please select a file to process.')
            continue
        
        else:
            # Use their chosen thread value, else default to 4
            if values['-THREADS-'] != '':
                chosen_threads = values['-THREADS-']
            else: 
                chosen_threads = 4
        
            # Run Taxomatcher on their chosen database
            chosen_csv = values['-IN_CSV-']

            if values['gbif']:
                gbif_main(input_csv = chosen_csv, threads = chosen_threads)
            elif values['ncbi']:
                if values['-EMAIL-'] == '':
                    sg.popup('NCBI requires entry of a user email to query.')
                    continue
                else: 
                    entrez_main(input_csv = chosen_csv, user_email=values['-EMAIL-'])

            # Lock in GBIF/NCBI choice
            window['gbif'].update(disabled=True)
            window['ncbi'].update(disabled=True)

            # Prompt the user to select a file destination
            window['-OUTPUT_MSG-'].update(visible=True)
            window['-DEST-'].update(visible=True)
            window['-DEST_BROWSE-'].update(visible=True)
   
    # Reveal save button if output path/name has been selected
    if event == '-DEST-' and values['-DEST-'] != '':
        window['Save'].update(visible=True)
    
    if event == 'Save':
        dest_dr = values['-DEST-']
        run_name = chosen_csv.split("/")[-1].split(".")[0]

        # Move output file from output directory to destination
        if values['gbif']:
            output_file = f"../output/{run_name}_gbif_output.tsv"
            final_location = os.path.join(dest_dr, os.path.basename(output_file))

            # Replace the files if they already exist in output
            if os.path.isfile(final_location):
                os.replace(final_location, output_file)

            shutil.move(output_file, dest_dr)
            
        elif values['ncbi']:
            output_file = f"../output/{run_name}_ncbi_output.tsv"
            fail_file = f"../output/{run_name}_ncbi_fails.tsv"
            dest_dr = values['-DEST-']

            # Move the files to the destination folder
            shutil.move(output_file, dest_dr)
            shutil.move(fail_file, dest_dr)
        
        sg.popup(f'Taxomatcher output has been saved to {dest_dr}')

        break


window.close()
import PySimpleGUI as sg
import os
import shutil
from gbif_matcher import main as gbif_main
from gbif_matcher import get_file_name as gbif_output_name

sg.theme('DarkAmber')   # set GUI color

# Define all the stuff inside GUI window
layout = [[sg.Text('Enter CSV file:')],
          
          [sg.Input(enable_events=True, key='-IN-'), sg.FileBrowse()], # Let GUI detect inputs and add browse button

          # [sg.Text('', size=(1))], # Add a line space
          [sg.Button('Run Taxomatcher')],
    
          [sg.Text('', size=(1))], # Add a line space
          [sg.Text('Select a destination folder:', visible=False, key='-OUTPUT_MSG-')],
          [sg.Input(enable_events=True, visible=False, key='-DEST-'), sg.FolderBrowse(visible=False, key='-DEST_BROWSE-')],

          [sg.Button('Save', visible=False), sg.Button('Exit')],
          ]

# Create the GUI window
window = sg.Window('Taxomatcher', layout, size = (450,220))

# Event Loop to process GUI events
while True:
    event, values = window.read()
    if event == sg.WIN_CLOSED or event == 'Exit': # if user closes window or clicks exit
        break

    # If they click run
    elif event == 'Run Taxomatcher':
        # Check if a file has been selected:
        if not values['-IN-']:
            sg.popup('Please select a file to process.')
            continue
        
        else:
            # Run Taxomatcher on selected file
            gbif_main(input_csv = values['-IN-'], threads = 4)

            # Prompt the user to select a file destination
            window['-OUTPUT_MSG-'].update(visible=True)
            window['-DEST-'].update(visible=True)
            window['-DEST_BROWSE-'].update(visible=True)
   
    # Reveal save button if output path/name has been selected
    if event == '-DEST-' and values['-DEST-'] != '':
        window['Save'].update(visible=True)
    
    if event == 'Save':
        dest_dr = values['-DEST-']

        output_file = gbif_output_name()

        # Move output file from output directory to destination
        shutil.move(output_file, dest_dr)
        
        sg.popup(f'Taxomatcher output has been saved to {dest_dr}')

        break


window.close()
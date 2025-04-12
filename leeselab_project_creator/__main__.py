import FreeSimpleGUI as sg
from leeselab_project_creator import generate_plate_layout, generate_worklist
from pathlib import Path


# main function to call the GUI
def main():
    # define the layout of the GUI
    layout = [
        [
            sg.Text(
                "This is the Leeselab laboratory planing assistant.\nPlease fill the fields to get started."
            )
        ],
        [sg.Text("Path to the list of samples:")],
        [
            sg.InputText(
                size=(40, 1),
                do_not_clear=True,
                key="input_file",
                tooltip="Please select a list of samples in Excel format.",
            ),
            sg.FileBrowse(
                file_types=(("Excel", ".xlsx"),),
                tooltip="Please select a list of samples in Excel format.",
            ),
        ],
        [sg.Text("Excel sheet name:")],
        [
            sg.InputText(
                default_text="Tabelle1",
                key="sheet_name",
                size=(40, 1),
                tooltip="Please write the name of the correct sheet here.",
            )
        ],
        [sg.Text("First sample can be found in:")],
        [
            [
                sg.Text("Column"),
                sg.InputText(
                    size=(10, 1),
                    tooltip="Indicate the column of the first sample (e.g. A)",
                    key="column",
                ),
                sg.Text("Row"),
                sg.InputText(
                    size=(10, 1),
                    tooltip="Indicate the row of the first samples (e.g. 1)",
                    key="row",
                ),
            ]
        ],
        [sg.Text("Number of extraction replicates:")],
        [
            sg.Spin(
                values=list(range(1, 501)),
                initial_value=1,
                size=(5, 1),
                key="extraction_replicates",
                tooltip="How many extraction replicates do you need?",
            )
        ],
        [sg.Text("Number of PCR replicates:")],
        [
            sg.Spin(
                values=list(range(1, 501)),
                initial_value=1,
                size=(5, 1),
                key="pcr_replicates",
                tooltip="How many PCR replicates do you need?",
            )
        ],
        [sg.Checkbox(text="Pool PCR replicates", default=False, key="pool_replicates")],
        [sg.Text("Primers to use:")],
        [
            sg.InputText(
                size=(40, 1),
                tooltip="Please add all primers you want to use seperated by comma (e.g. fwh2, rbcl)",
                key="markers",
            )
        ],
        [sg.Text("Available 1st step tags:")],
        [
            sg.InputText(
                size=(40, 1),
                tooltip="Please add the number of the available first step tags (e.g. 1, 2, 5, 6)",
                key="primers",
            )
        ],
        [sg.Checkbox(text="Suggest pooling strategy", default=True, key="pooling")],
        [sg.Text("Project name")],
        [
            sg.InputText(
                size=(40, 1),
                tooltip="Please enter a project name here",
                key="project_name",
            )
        ],
        [sg.Text("Select an output folder")],
        [
            sg.InputText(
                size=(40, 1),
                tooltip="Please select an output folder to save the generated sheets",
                key="output_folder",
            ),
            sg.FolderBrowse(
                tooltip="Please select an output folder to save the generated sheets"
            ),
        ],
        [
            sg.Button(
                "Generate sample list",
                button_color="red",
                size=(42, 2),
                disabled=True,
                key="generate",
            )
        ],
        [
            sg.Multiline(
                size=(45, 5),
                key="output_text",
                autoscroll=True,
                default_text="Please fill all fields to generate a project.\n",
            )
        ],
        [
            sg.Button(
                "Cancel",
            )
        ],
    ]

    # define the window to be called
    window = sg.Window("Leeselab project creator", layout)
    # was the output already updated?
    updated_out = False
    # run the main loop
    while True:
        # read the windows values and events every 100 ms
        event, values = window.read(timeout=200)

        # remove those values since they are not needed, easier to check to update the button
        try:
            del values["Browse"]
            del values["Browse0"]
        except TypeError:
            pass

        # if the close button is pressed, break the loop, close the window
        if event == sg.WIN_CLOSED or event == "Cancel":
            break

        # find out if we are ready for processing
        ready = [i for i in values.values()]
        ready = ready[:6] + ready[7:8] + ready[10:]

        # enable generate button when all fields are filled
        if all(ready):
            window["generate"].update(disabled=False)
            if not updated_out:
                window["output_text"].print("Ready to generate project files.")
                updated_out = True
            else:
                updated_out = True
        else:
            window["generate"].update(disabled=True)
            updated_out = False

        if event == "generate":
            # generate the plate layout
            generate_plate_layout.generate_plate_layout(
                input_file=Path(values["input_file"]),
                sheet=values["sheet_name"],
                col=values["column"],
                row=values["row"],
                extraction_replicates=values["extraction_replicates"],
                output_folder=Path(values["output_folder"]),
                project_name=values["project_name"],
            )

            # generate the worklist
            generate_worklist.generate_worklist(
                output_path=Path(values["output_folder"]),
                project=values["project_name"],
                available_primers=values["primers"],
                pcr_replicates=values["pcr_replicates"],
                markers=values["markers"],
                pool_pcr_replicates=values["pool_replicates"],
            )
            window["output_text"].print(
                "Plate layout and worklist saved to output folder."
            )

    window.Close()


# run only if called as a toplevel script
if __name__ == "__main__":
    main()

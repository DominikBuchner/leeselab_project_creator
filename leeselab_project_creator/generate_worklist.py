import pandas as pd
import math


def generate_worklist(output_path, project, available_primers, pcr_replicates, markers):
    # generate the correct path to the file
    input_file = "{}_plate_layout.xlsx".format(project)
    input_file = output_path.joinpath(input_file)

    # read the data
    data = pd.read_excel(input_file, sheet_name="plate_layout", skiprows=1)

    # collect all worklists here for the final output
    worklists = []

    # collect the plate letters
    plates = list(data["Lysis,\nExtraction,\nPCR plate"].unique())
    general_worklist = pd.DataFrame(
        columns=[
            "plate",
            "aliquoted",
            "lysis",
            "inhibitor removal",
            "lysate distribution",
            "extraction",
            "gel extraction",
        ]
    )

    # general worklist is done
    general_worklist["plate"] = plates

    # add the library x pcr replicates worklist
    # calculate the number of librarys needed
    available_primers = available_primers.split(",")
    markers = markers.split(",")
    number_of_extraction_plates = len(plates)
    number_of_pcr_plates = number_of_extraction_plates * pcr_replicates
    number_of_librarys = math.ceil(number_of_pcr_plates / len(available_primers))
    plates_per_library = int(number_of_pcr_plates / number_of_librarys)
    optimal_primer_order = [
        1,
        5,
        9,
        13,
        17,
        21,
        2,
        6,
        10,
        14,
        18,
        22,
        3,
        7,
        11,
        15,
        19,
        23,
        4,
        8,
        12,
        16,
        20,
        24,
    ]

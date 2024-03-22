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

    # calculate everything marker wise and concat the marker dfs in the end

    # add the library x pcr replicates worklist
    # calculate the number of librarys needed
    available_primers = available_primers.split(",")
    markers = markers.split(",")
    number_of_extraction_plates = len(plates)
    number_of_pcr_plates = number_of_extraction_plates * pcr_replicates * len(markers)
    librarys_per_marker = math.ceil(
        (number_of_pcr_plates / len(markers)) / len(available_primers)
    )
    plates_per_library_per_marker = int(
        (number_of_pcr_plates / librarys_per_marker) / len(markers)
    )

    print(
        number_of_extraction_plates,
        number_of_pcr_plates,
        librarys_per_marker,
        plates_per_library_per_marker,
    )

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

    # generate the output dataframe
    worklist_dataframe = pd.DataFrame()
    # add the source plate
    worklist_dataframe["source_plate"] = [
        plate_letter for plate_letter in plates for _ in range(pcr_replicates)
    ] * len(markers)

    # add the respective librarys
    worklist_dataframe["library"] = [
        library
        for library in range(1, librarys_per_marker * len(markers) + 1)
        for _ in range(plates_per_library_per_marker)
    ]

    # calculate the optimal primer order for the 1st pcr
    print(worklist_dataframe)

    #### code still broken for uneven number of plates in the last library. Have to cut somehow. #####

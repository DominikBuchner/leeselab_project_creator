import pandas as pd
import math


# funtion to calculate one worklist per marker
def worklist_per_marker(
    available_primers,
    extraction_plates,
    pcr_replicates,
    primer_name,
    optimal_primers,
    starting_library_number,
):
    # generate the output dataframe
    output_worklist = pd.DataFrame()
    output_worklist["source plate"] = [
        extraction_plate
        for extraction_plate in extraction_plates
        for _ in range(pcr_replicates)
    ]

    # find the number of pcr plates for this marker
    pcr_plate_count = len(output_worklist)

    # find the number of librarys to distribute primers evenly
    available_primers = available_primers.split(",")
    available_primers = [int(primer) for primer in available_primers]
    number_of_primers = len(available_primers)

    while True:
        if pcr_plate_count <= number_of_primers:
            library_count = 1
            break
        elif pcr_plate_count % number_of_primers == 0:
            library_count = int(pcr_plate_count / number_of_primers)
            break
        else:
            number_of_primers -= 1

    # calculate the number of plates per library
    plates_per_library = int(pcr_plate_count / library_count)

    # generate the library column
    library_column = []

    for i in range(library_count):
        for _ in range(plates_per_library):
            library_column.append(starting_library_number)
        starting_library_number += 1

    output_worklist["library"] = library_column

    primers = []

    # find the optimal primers for maximum library diversity
    for primer in optimal_primers:
        if len(primers) < plates_per_library:
            if primer not in primers and primer in available_primers:
                primers.append(primer)
        else:
            break

    primers = primers * library_count
    primers = [
        "{} - {}".format(primer_name, primer_number) for primer_number in primers
    ]

    output_worklist["tagging primer"] = primers
    output_worklist["1st pcr"] = ""
    output_worklist["clean up"] = ""
    output_worklist["2nd pcr"] = ""
    output_worklist["normalization"] = ""
    output_worklist["pooling"] = ""

    return output_worklist, starting_library_number


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

    # optimal primer order for maximizing library diversity
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

    # gather the worklists here
    worklists = []

    # calculate everything marker wise and concat the marker dfs in the end
    next_library = 1
    markers = markers.split(",")

    for primer in markers:
        worklist, next_library = worklist_per_marker(
            available_primers,
            plates,
            pcr_replicates,
            primer,
            optimal_primer_order,
            next_library,
        )

        worklists.append(worklist)

    # concat the individual worklist to have a working table
    working_table = pd.concat(worklists, axis=0)
    print(working_table)

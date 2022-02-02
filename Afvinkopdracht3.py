from translation import code


def main():
    dna_1 = "AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGA" \
            "GTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG"
    dna_2 = complimentary(dna_1)
    test_sequence(dna_1)
    test_sequence(dna_2)


def translation(sequence, reading_frame):
    """Checks every 3 characters of string, if matched in dict
    adds value to protein string.

    :param sequence -- string, DNA provided in program.
    :param reading_frame -- int, user input.
    :return protein -- string, returns a formatted protein sequence"""
    sequence = sequence.lower()
    length = len(sequence) // 3
    protein = ""
    for dna_index in range(length):
        start_index = (dna_index * 3) + reading_frame
        end_index = start_index + 3
        dna_segment = sequence[start_index:end_index]
        if dna_segment in code:
            protein += code[dna_segment]
    return protein


def complimentary(sequence):
    """Makes the provided string of DNA complimentary and reversed."""
    return sequence.replace("A", "t") \
        .replace("T", "a")\
        .replace("C", "g")\
        .replace("G", "c").upper()[::-1]


def test_sequence(sequence):
    """"Loops through  frames for different mRNA sequences.
    :param sequence -- string, DNA provided in program.
    :returns print statement"""
    for frame in range(3):
        protein = translation(sequence, frame)
        print(f"The protein for frameshift {frame + 1} is: \n"
              f"{protein}\n")


if __name__ == '__main__':
    main()

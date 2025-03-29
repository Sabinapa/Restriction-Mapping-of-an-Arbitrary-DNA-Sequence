def read_dna_file(filepath):
    """
    Reads a DNA sequence from a text file and returns it as a string.

    :param filepath: Path to the DNA text file.
    :return: DNA sequence as a string.
    """

    with open(filepath, 'r') as file:
        dna_sequence = file.read().strip()
    return dna_sequence


def find_restriction_sites(dna_sequence, restriction_sequence):
    """
    Finds all occurrences of the restriction sequence within the DNA sequence.

    :param dna_sequence: The DNA sequence to search in.
    :param restriction_sequence: The restriction sequence to find.
    :return: A list of indices where the restriction sites start.
    """
    sites = []
    index = dna_sequence.find(restriction_sequence)
    while index != -1:
        sites.append(index)
        index = dna_sequence.find(restriction_sequence, index + 1)
    return sites


# Example usage
if __name__ == "__main__":
    dna1 = read_dna_file('DNK1.txt')
    dna2 = read_dna_file('DNK2.txt')
    dna3 = read_dna_file('DNK3.txt')

    #print("DNA1 sequence:", dna1)
    #print("DNA2 sequence:", dna2)
    #print("DNA3 sequence:", dna3)

    restriction_seq = "GAATTC"  # primer restrikcijskega reza (lahko spremeniš)

    sites = find_restriction_sites(dna2, restriction_seq)
    print(f"Restriction sequence '{restriction_seq}' found at indices:", sites)

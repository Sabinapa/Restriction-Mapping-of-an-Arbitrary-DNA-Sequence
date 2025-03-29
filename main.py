def read_dna_file(filepath):
    """
    Reads a DNA sequence from a text file and returns it as a string.

    :param filepath: Path to the DNA text file.
    :return: DNA sequence as a string.
    """

    with open(filepath, 'r') as file:
        dna_sequence = file.read().strip()
    return dna_sequence


# Example usage
if __name__ == "__main__":
    dna1 = read_dna_file('DNK1.txt')
    dna2 = read_dna_file('DNK2.txt')
    dna3 = read_dna_file('DNK3.txt')

    print("DNA1 sequence:", dna1)
    #print("DNA2 sequence:", dna2)
    #print("DNA3 sequence:", dna3)

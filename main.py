from itertools import combinations

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


def calculate_multiset_of_distances(restriction_sites):
    """
    Calculates multiset of distances between all pairs of restriction sites.

    :param restriction_sites: List of indices of restriction sites.
    :return: Sorted list representing the multiset of distances.
    """
    distances = []
    for site1, site2 in combinations(restriction_sites, 2):
        distance = abs(site2 - site1)
        distances.append(distance)
    distances.sort()
    return distances

def get_pairwise_distances(points):
    """Return sorted multiset of pairwise distances."""
    distances = sorted(abs(b - a) for a, b in combinations(points, 2))
    return distances

def naive_restriction_sites_reconstruction(distance_multiset):
    width = max(distance_multiset)
    potential_sites = range(1, width)
    solution = None

    for subset_size in range(1, len(potential_sites)+1):
        for subset in combinations(potential_sites, subset_size):
            candidate = [0] + list(subset) + [width]
            if get_pairwise_distances(candidate) == sorted(distance_multiset):
                solution = candidate
                return sorted(solution)

    return "No solution found"

def partial_digest(L):
    width = max(L)
    L.remove(width)
    X = [0, width]
    place(L, X, width)

def place(L, X, width):
    if not L:
        X.sort()
        print("Solution found:", X)
        return

    y = max(L)
    Dy = [abs(y - x) for x in X]

    # Try placing y
    if all(Dy.count(dist) <= L.count(dist) for dist in Dy):
        for dist in Dy:
            L.remove(dist)
        X.append(y)
        place(L, X, width)
        X.remove(y)
        for dist in Dy:
            L.append(dist)

    # Try placing width - y
    Dy = [abs(width - y - x) for x in X]
    if all(Dy.count(dist) <= L.count(dist) for dist in Dy):
        for dist in Dy:
            L.remove(dist)
        X.append(width - y)
        place(L, X, width)
        X.remove(width - y)
        for dist in Dy:
            L.append(dist)


# Example usage
if __name__ == "__main__":
    dna1 = read_dna_file('DNK1.txt')
    dna2 = read_dna_file('DNK2.txt')
    dna3 = read_dna_file('DNK3.txt')

    #print("DNA1 sequence:", dna1)
    #print("DNA2 sequence:", dna2)
    #print("DNA3 sequence:", dna3)

    restriction_seq = "GTGTG"  # primer restrikcijskega reza (lahko spremeniš)

    sites = find_restriction_sites(dna1, restriction_seq)
    print(f"Restriction sequence '{restriction_seq}' found at indices:", sites)

    multiset_distances = calculate_multiset_of_distances(sites)
    print("Multiset of distances between restriction sites:", multiset_distances)

    reconstruction = naive_restriction_sites_reconstruction(multiset_distances)
    print("Reconstructed restriction sites:", reconstruction)

    print("Starting Branch and Bound approach:")
    partial_digest(multiset_distances.copy())

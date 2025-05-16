import random

def generate_random_vcf(filename='sample.vcf', total_variants=1000):
    # Chromosomes 1-22 plus X and Y
    chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']

    # Randomly assign variant counts per chromosome that sum to total_variants
    # Method: generate 25 random weights, normalize to total_variants
    weights = [random.random() for _ in chromosomes]
    total_weight = sum(weights)
    counts = [int(w / total_weight * total_variants) for w in weights]

    # Adjust counts to sum exactly to total_variants
    diff = total_variants - sum(counts)
    for i in range(abs(diff)):
        index = i % len(counts)
        counts[index] += 1 if diff > 0 else -1

    # Prepare VCF header lines
    header = [
        "##fileformat=VCFv4.2",
        "##source=RandomVariantGenerator",
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    ]

    # Generate variants with random positions and alleles
    with open(filename, 'w') as f:
        for line in header:
            f.write(line + '\n')

        # For each chromosome, generate the specified number of variants
        for chrom, count in zip(chromosomes, counts):
            for _ in range(count):
                pos = random.randint(1, 1_000_000)  # position 1 to 1,000,000
                var_id = '.'  # no ID
                ref = random.choice(['A', 'C', 'G', 'T'])
                alt = random.choice([b for b in ['A', 'C', 'G', 'T'] if b != ref])
                qual = round(random.uniform(10, 60), 2)
                filt = 'PASS'
                info = f"DP={random.randint(10, 100)}"
                line = f"{chrom}\t{pos}\t{var_id}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}"
                f.write(line + '\n')

    print(f"Generated {total_variants} variants in '{filename}'.")

if __name__ == "__main__":
    generate_random_vcf()


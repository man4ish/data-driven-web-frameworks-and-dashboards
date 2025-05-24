import random

def random_sequence(length):
    return ''.join(random.choices('ACGT', k=length))

def random_quality(length):
    # Generate dummy quality scores (Phred+33 ASCII range from '!' to 'I')
    return ''.join(random.choices([chr(q) for q in range(33, 73)], k=length))

def generate_fastq(filename, num_reads=10, read_length=100):
    with open(filename, 'w') as f:
        for i in range(1, num_reads + 1):
            seq = random_sequence(read_length)
            qual = random_quality(read_length)
            f.write(f"@read{i}\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write(f"{qual}\n")

if __name__ == "__main__":
    generate_fastq("sample_reads.fastq", num_reads=20, read_length=150)


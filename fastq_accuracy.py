import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

def phred_quality_to_accuracy(qual):
    """converts Phred quality values into base precision (1-P)"""
    return [1 - 10 ** (-q / 10) for q in qual]

def compute_accuracy_stats(fastq_file):
    """reads a FASTQ file and calculates the precision of the bases considering length differences"""
    accuracies_by_position = []

    # to iterate over each FASTQ file read
    for record in SeqIO.parse(fastq_file, "fastq"):
        acc = phred_quality_to_accuracy(record.letter_annotations["phred_quality"])

        # to ensure that the storage structure fits the length of the read
        if len(accuracies_by_position) < len(acc):
            accuracies_by_position.extend([[] for _ in range(len(acc) - len(accuracies_by_position))])

        # to store the precision of each base in its corresponding position
        for i, value in enumerate(acc):
            accuracies_by_position[i].append(value)

    max_length = len(accuracies_by_position)

    # to calculate mean precision and standard deviation per position, ignoring empty positions
    mean_per_site = np.array([np.mean(pos) if pos else np.nan for pos in accuracies_by_position])
    std_per_site = np.array([np.std(pos) if pos else np.nan for pos in accuracies_by_position])

    # to calculate the overall mean and standard deviation of precision
    mean_accuracy = np.nanmean(mean_per_site) # I use the 'nan' to ignore empty data
    std_accuracy = np.nanstd(mean_per_site)

    print(f"Mean base call accuracy: {mean_accuracy}")
    print(f"Standard deviation of base call accuracy: {std_accuracy}")

    return mean_per_site, std_per_site, max_length


def plot_accuracy(mean_per_site, std_per_site, max_length, n):
    """generates a graph of the average precision by position with 95% and 99% confidence intervals"""
    ci_95 = 1.96 * (std_per_site / np.sqrt(n))  # 95% confidence interval
    ci_99 = 2.58 * (std_per_site / np.sqrt(n))  # 99% confidence interval

    plt.figure(figsize=(10, 5))
    plt.plot(mean_per_site, label='Mean Accuracy')
    plt.fill_between(range(max_length), mean_per_site - ci_95, mean_per_site + ci_95, color='blue', alpha=0.2,
                     label='95% CI')
    plt.fill_between(range(max_length), mean_per_site - ci_99, mean_per_site + ci_99, color='red', alpha=0.2,
                     label='99% CI')

    plt.xlabel("Base Position")
    plt.ylabel("Accuracy")
    plt.title("Mean Base Call Accuracy Per Site")
    plt.legend()
    plt.show()


def main():
    """it manages the arguments and execution of base precision analysis"""
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq", type=str, help="Path to the FASTQ file")
    args = parser.parse_args()

    # to validate that the specified FASTQ file exists and is not a directory
    if not os.path.isfile(args.fastq):
        print("the specified FASTQ file does not exist or is a directory")
        return

    # to calculate accuracy statistics and graph results
    mean_per_site, std_per_site, max_length = compute_accuracy_stats(args.fastq)
    plot_accuracy(mean_per_site, std_per_site, max_length, max_length)


if __name__ == "__main__":
    main()

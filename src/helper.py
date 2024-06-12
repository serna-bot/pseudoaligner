import matplotlib
from matplotlib import pyplot as plt
from tqdm import tqdm
from Bio import SeqIO
import os.path
from pathlib import Path


def create_histogram(fig, x, y, xlabel, ylabel, title, num_bars = 1, legend = None):
    plt.figure(fig)
    if num_bars > 1:
        for idx, y1 in enumerate(y):
            plt.bar(x, y1, label = legend[idx])
        plt.legend() 
    else:
        plt.bar(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

def create_plot(fignum, x, y, xlabel, ylabel, title, xticks = None):
    plt.figure(fignum)
    plt.plot(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if xticks:
        plt.xticks(x, xticks)

def get_ground_truth(src):
    equivalence_classes = {}
    for record in tqdm(SeqIO.parse(src, 'fasta')):
        identifier = record.id
        equiv_class = identifier[identifier.find('/') + 1:identifier.find(';')]
        if equiv_class in equivalence_classes:
            equivalence_classes[equiv_class] += 1
        else:
            equivalence_classes[equiv_class] = 1
    file_path = "/Users/serenaong/Documents/csc121_project/ground_truth.txt"
    file = Path(file_path)
    if not file.is_file():
        f = open(file_path, "a")
        for key, value in equivalence_classes.items():
            f.write(f"{key}:{value}\n")
        f.close()
    return equivalence_classes

def save_to_file(equivalence_classes, kmer_len, skipping):
    file_name = "/Users/serenaong/Documents/csc121_project/equiv_classes" + str(kmer_len)
    if skipping:
        file_name += "skipping"
    file_output = ""
    f = open(file_name, "w")
    for key, value in equivalence_classes.items():
        file_output += f"{key}:{value}\n"
    f.write(file_output)
    f.close()
    print("Saved to file.")

def get_reads_stats(src):
    f = open(src, "r")
    lines = f.readlines()
    max_equiv_class_len = 0
    num_equiv_classes_dict = {}
    num_in_na, percent_na = 0, 0.0
    total_reads = 0
    for line in lines:
        sep = line.find(":")
        equiv_classes = line[:sep]
        num_reads = int(line[sep + 1:])
        total_reads += num_reads
        if equiv_classes != "NA":
            num_equiv_classes = equiv_classes.count(',') + 1
            max_equiv_class_len = max(max_equiv_class_len, num_equiv_classes)
            if num_equiv_classes in num_equiv_classes_dict:
                num_equiv_classes_dict[num_equiv_classes]['count'] += 1
                num_equiv_classes_dict[num_equiv_classes]['reads'] += num_reads
            else:
                num_equiv_classes_dict[num_equiv_classes] = {'count': 1, 'reads': num_reads}
        else:
            num_in_na = num_reads
    percent_na = num_in_na / total_reads
    return {'num_in_na': num_in_na, 'percent_na': percent_na, 'num_equiv_classes_dict': num_equiv_classes_dict, 'max_equiv_class_len': max_equiv_class_len}
    
import matplotlib
from matplotlib import pyplot as plt
from tqdm import tqdm
from Bio import SeqIO
import os.path
from pathlib import Path



def create_figure(equiv_classes):
    pass

def create_plot(plots, equiv_classes, single_figure=False):
    #check if the fig file is present, else create one
    pass

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

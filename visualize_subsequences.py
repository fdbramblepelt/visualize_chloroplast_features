#!/usr/bin/env python3
import csv
import argparse
import pandas as pd
import random
from collections import Counter
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

TEXT_SIZE = 0

def split_header(header):
    try:
        header1, hairpin_id, sequence_id, location, remaining = header.split('-', 4)
        genus, species, region = header1.split('_')
        # Use regular expressions to split based on the first occurrence of ':'
        neighbor_gene_match = re.match(r"(.+):(.+)-(.+):(.+)-(.+)", remaining.split('_')[0])

        if not neighbor_gene_match:
            raise ValueError("Gene information is not in the expected format.")

        neighbor_gene_name, neighbor_gene_distance, second_neighbor_gene_name, second_neighbor_gene_distance, tribe = neighbor_gene_match.groups()
        subfamily, major_clade, photosynthetic_pathway = remaining.split('_')[1:]

        neighbor_gene = {
            "name": neighbor_gene_name,
            "distance": neighbor_gene_distance
        }
        second_neighbor_gene = {
            "name": second_neighbor_gene_name,
            "distance": second_neighbor_gene_distance
        }
    except ValueError:
        print(f"Error: Header '{header}' does not follow the expected format.")
        return None
    return {
        "sequence_id": sequence_id,
        "genus": genus,
        "species": species,
        "hairpin_id": region + "_" +hairpin_id,
        "region": region,
        "local_hairpin_id": hairpin_id,
        "start_pos": location.split(':')[0],
        "end_pos": location.split(':')[1],
        "length": int(location.split(':')[1]) - int(location.split(':')[0]),    
        "neighbor_gene_name": neighbor_gene["name"],
        "neighbor_gene_distance": neighbor_gene["distance"],
        "second_neighbor_gene_name": second_neighbor_gene["name"],
        "second_neighbor_gene_distance": second_neighbor_gene["distance"],
        "tribe": tribe,
        "subfamily": subfamily,
        "major_clade": major_clade,
        "photosynthetic_pathway": photosynthetic_pathway
    }

def parse_data(rows):
    parsed_data = []

    for row in rows:
        columns = row.split(',')
        header_info = split_header(columns[0])
        if header_info is None:
            continue

        alias_count = int(columns[1])
        group_number = columns[2]
        scores = {f"G{i+1}": int(columns[i+3]) for i in range(len(columns) - 3)}

        parsed_data.append({
            "header_info": header_info,
            "alias_count": alias_count,
            "group_number": group_number,
            "scores": scores
        })

    return parsed_data

################################################################################

def extract_score_keys(parsed_data):
    """Extract and sort score keys from parsed data, replacing scores with 0 when group_number{i} == group_number{j}."""
    score_keys = set()
    for data in parsed_data:
        for key in data["scores"].keys():
            score_keys.add(key)
            if data["group_number"] == key:
                data["scores"][key] = 10
    return sorted(score_keys, key=lambda x: int(x[1:]))

def extract_labels_and_groups(parsed_data):
    """Extract genus-species labels and group numbers for sorting."""
    genus_species_labels = sorted(
        (
            int(data['group_number'][1:]), 
                f"{data['header_info']['genus']}_{data['header_info']['species']} | {data['header_info']['sequence_id']} | {data['header_info']['hairpin_id']} | "
                f"len: {data['header_info']['length']} | {data['group_number']}\n2nd-neighbor: {data['header_info']['second_neighbor_gene_name']} #-in-db: {data['alias_count']} | "
                f"{data['header_info']['tribe']}, {data['header_info']['subfamily']}, {data['header_info']['major_clade']} | {data['header_info']['photosynthetic_pathway']}\n"
        )
        for data in parsed_data
    )
    labels = [label for _, label in genus_species_labels]
    group_numbers = [data['group_number'] for data in parsed_data]
    return labels, group_numbers

def build_score_matrix(parsed_data, score_keys, genus_species_labels, group_numbers):
    """Build a score matrix from parsed data."""
    score_matrix = np.zeros((len(genus_species_labels), len(score_keys)))
    for i, genus_species in enumerate(genus_species_labels):
        group_number = f"G{i+1}"
        for data in parsed_data:
            if data['group_number'] == group_number:
                for j, score_key in enumerate(score_keys):
                    score_matrix[i, j] = data["scores"].get(score_key, 0)
    return score_matrix

################################################################################

def calculate_text_size(score_keys, genus_species_labels):
    # Example logic to calculate text size based on the number of score keys and genus-species labels
    num_score_keys = len(score_keys)
    num_genus_species = len(genus_species_labels)

    if max(num_score_keys, num_genus_species) > 200:
        text_size = 2
    elif max(num_score_keys, num_genus_species) > 100:
        text_size = 4
    elif max(num_score_keys, num_genus_species) > 50:
        text_size = 6
    elif max(num_score_keys, num_genus_species) > 25:
        text_size = 8
    else:
        text_size = 12

    print(f"Debug: Calculated text size: {text_size}")
    return text_size

################################################################################

def plot_heatmap(score_matrix, score_keys, genus_species_labels, group_numbers, output_file_path):
    """Plot a heatmap of the score matrix."""
    global TEXT_SIZE
    #TEXT_SIZE = calculate_text_size(score_matrix.shape[0], score_matrix.shape[1])
    cmap = plt.cm.viridis
    norm = mcolors.Normalize(vmin=score_matrix.min(), vmax=score_matrix.max())

    fig, ax = plt.subplots(figsize=(15, 10))
    cax = ax.matshow(score_matrix, cmap=cmap, norm=norm)
    add_color_bar(fig, cax)
    annotate_cells(ax, score_matrix, group_numbers, norm)
    set_axis_labels(ax, score_keys, genus_species_labels)
    set_title(ax, output_file_path)
    save_and_show_plot(fig, output_file_path)

def visualize_scores(parsed_data, output_file_path):
    score_keys = extract_score_keys(parsed_data)
    genus_species_labels, group_numbers = extract_labels_and_groups(parsed_data)
    score_matrix = build_score_matrix(parsed_data, score_keys, genus_species_labels, group_numbers)
    plot_heatmap(score_matrix, score_keys, genus_species_labels, group_numbers, output_file_path)

################################################################################

def annotate_cells(ax, score_matrix, group_numbers, norm):
    """Annotate each cell in the heatmap with the score."""
    global TEXT_SIZE
    fontsize = TEXT_SIZE
    alpha = 0.7
    text = ""
    for i in range(score_matrix.shape[0]):
        max_score = max(score_matrix[i, j] for j in range(score_matrix.shape[1]) if group_numbers[i] != group_numbers[j])
        for j in range(score_matrix.shape[1]):
            score = score_matrix[i, j]
            # Determine text color based on score comparison
            if group_numbers[i] == group_numbers[j]:
                color = 'blue'
                fontweight = 'bold'
                text = "NA"
            elif score >= max_score:
                #print("Score is big :)") DEBUG - makes sure max scores are being registered from data.
                color = 'magenta'
                fontweight = 'bold'
                text = int(score)
            else:
                color = plt.cm.viridis(norm(score * 0.4))  # Darker shade for other scores
                fontweight = 'normal'
                text = int(score)
            
            ax.text(j, i, f'{text}', va='center', ha='center', color=color, fontsize=fontsize, fontweight=fontweight, alpha=alpha)

################################################################################

def add_color_bar(fig, cax):
    """Add a color bar to the plot."""
    cbar = fig.colorbar(cax, orientation='horizontal', pad=0.1)
    cbar.ax.tick_params(labelsize=8)

def set_axis_labels(ax, score_keys, genus_species_labels):
    """Set labels for the axes."""
    fontsize=TEXT_SIZE
    ax.set_xticks(np.arange(len(score_keys)))
    ax.set_yticks(np.arange(len(genus_species_labels)))
    ax.set_xticklabels([f"G{i+1}" for i in range(len(score_keys))], rotation=45, ha='right', fontsize=((int(fontsize) / 3) * 2))
    ax.set_yticklabels(genus_species_labels, fontsize=((int(fontsize) / 3) * 2))
    ax.set_xlabel('Score Keys')
    ax.set_ylabel('Genus - Species')

def set_title(ax, output_file_path):
    """Set the title of the plot based on the output file path."""
    global TEXT_SIZE
    match = re.search(r".+subsequence_scores\.(.+)\.fasta.+", output_file_path)
    if match:
        extracted_phrase = match.group(1)
        ax.set_title(f'UNIQUE Sequence Substring Similarity Scores by Genus - Species (neighbor: {extracted_phrase})', fontsize=TEXT_SIZE)
    else:
        ax.set_title('Heatmap of Substring Sequence Similarity Scores by Genus - Species', fontsize=TEXT_SIZE)

def save_and_show_plot(fig, output_file_path):
    """Save the plot to a file and display it."""
    fig.set_size_inches(20, 20)
    plt.tight_layout()
    plt.savefig(f"{output_file_path}.heatmap.png", dpi=300)  # Increase the dpi for higher quality
    plt.show()

################################################################################    

def main():
    parser = argparse.ArgumentParser(description='Process sequence data and aliases.')
    parser.add_argument('-i', '--input', default="default_csv", help='input CSV file')
    parser.add_argument('-od', '--output_dir', default='default_output_directory', help='output directory for files')
    args = parser.parse_args()

    if args.input == "default_csv":
        args.input = "/Users/fdbramblepelt/Desktop/Lab/Hairpin_Chloroplast/feb17_nearest_neighbor_data/pretty_subsequence_scores/pretty.subsequence_scores.atpE.fasta.csv.unique_sequences.csv"
    input_dir = '/'.join(args.input.split('/')[:-1])
    if args.output_dir == "default_output_directory":
        args.output_dir = input_dir + "/heatmaps"

    output_file_name = "HEATMAP." + args.input.split('/')[-1]
    output_file_path = f"{args.output_dir}/{output_file_name}"

    print(f"Debug: Reading input file {args.input}")
    print(f"Debug: Writing output file {output_file_path}")
    
    with open(args.input, 'r') as file:
        data = file.read()

    # Split the data into rows
    rows = data.split('\n')
    
    # Print the first 5 rows
    #for row in rows[:5]:
    #    print(row)

    parsed_data = parse_data(rows)
    global TEXT_SIZE
    TEXT_SIZE = calculate_text_size(extract_score_keys(parsed_data), extract_labels_and_groups(parsed_data)[0])
    visualize_scores(parsed_data, output_file_path)
    print(f"Done! Saved to {output_file_path}.heatmap.png")

if __name__ == "__main__":
    main()
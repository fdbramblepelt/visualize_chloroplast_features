#!/usr/bin/env python3

# Required packages:
# - matplotlib
# - bcbio-gff
# Install using: pip install matplotlib bcbio-gff

# Import necessary libraries

import random
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
import pandas as pd
import csv
import argparse

def parse_gff3_file(file_path, sequence_id):
    """
    Parse the GFF3 file to extract gene features for a specific sequence ID.
    """
    features = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue  # Skip comment lines

            parts = line.strip().split(',')
            if len(parts) != 9:
                continue  # Skip malformed lines

            seq_id, source, feature_type, start, end, score, strand, phase, attributes = parts

            # Check if the line corresponds to the requested sequence ID
            if not seq_id.lower().startswith(sequence_id.lower()):
                continue  # Skip lines that don't match the sequence ID

            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                key, value = attr.split('=')
                attr_dict[key] = value

            feature = {
                "type": feature_type,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "attributes": attr_dict
            }
            features.append(feature)
            print(f"Detected feature: {feature}")  # Debug statement

    print(f"Total features detected: {len(features)}")  # Debug statement
    return features




def parse_regions_file(region_file_path, sequence_id):
    """
    Parse the regions file to extract region sizes (LSC, SSC, IRs) for a specific sequence ID.
    """
    regions = {"LSC": None, "SSC": None, "IR": []}
    seen_ir = set()  # To track already seen IR regions

    with open(region_file_path, 'r') as file:
        for line in file:
            if line.startswith("#") or not line.strip():
                continue  # Skip comments and empty lines
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue  # Skip malformed lines
            
            # Check if the line corresponds to the requested sequence ID
            if not parts[0].startswith(sequence_id):
                continue  # Skip lines that don't match the sequence ID
            
            region_type = parts[2].lower()
            start = int(parts[3])
            end = int(parts[4])
            
            if region_type == "lsc":
                regions["LSC"] = (start, end)
            elif region_type == "ssc":
                regions["SSC"] = (start, end)
            elif region_type in ["ira", "irb"]:
                # Only add unique IR regions
                ir_key = (start, end)
                if ir_key not in seen_ir:
                    seen_ir.add(ir_key)
                    regions["IR"].append(ir_key)
    
    # Debug output
    print(f"Parsed regions for {sequence_id}: {regions}")
    
    return regions




def calculate_region_positions(regions):
    """
    Calculate the start and end positions for each region based on the parsed region lengths.
    """
    lsc_start, lsc_end = regions["LSC"]
    ir_length = regions["IR"][0][1] - regions["IR"][0][0]  # Assuming both IRa and IRb have the same length
    ssc_length = regions["SSC"][1] - regions["SSC"][0]

    ira_start = lsc_end + 1
    ira_end = ira_start + ir_length - 1
    ssc_start = ira_end + 1
    ssc_end = ssc_start + ssc_length - 1
    irb_start = ssc_end + 1
    irb_end = irb_start + ir_length - 1

    return {
        "LSC": (lsc_start, lsc_end),
        "IRa": (ira_start, ira_end),
        "SSC": (ssc_start, ssc_end),
        "IRb": (irb_start, irb_end)
    }





def organize_features(features, regions):
    """
    Organize features into genes and tRNAs.
    """
    genes = []
    tRNAs = []
    
    # Process genes and tRNAs as before
    for feature in features:
        if feature['type'] == "gene":
            print(f"Gene: {feature}")
            genes.append({
                "name": feature['attributes'].get("Name", "Unknown"),
                "start": feature['start'],
                "id": feature['attributes'].get("ID", "Unknown"),
                "end": feature['end'],
                "strand": feature['strand']
            })
        elif feature['type'] == "tRNA":
            print(f"tRNA: {feature}")
            tRNAs.append({
                "name": feature['attributes'].get("Name", "Unknown"),
                "id": feature['attributes'].get("ID", "Unknown"),
                "start": feature['start'],
                "end": feature['end'],
                "strand": feature['strand']
            })

    return regions, genes, tRNAs


####### COLOR GENES BY FUNCTION #######
def load_gene_colors():
    df = pd.read_csv("/Users/fdbramblepelt/Desktop/Lab/Hairpin_Chloroplast/Chloroplast_Vis/visualize_chloroplast_features/chloroplast_gene.csv")
    return {row['gene']: row['hex_code'].strip() for _, row in df.iterrows()}

def get_gene_color(gene_name, color_map):
    """
    Find the color for a gene based on partial name matching.
    """
    # Remove 'trn' prefix for tRNA matching to avoid double matches
    for ref_gene, color in color_map.items():
        if gene_name.startswith('trn') and ref_gene.startswith('trn'):
            if gene_name.split('-')[0] == ref_gene.split('-')[0]:
                return color
        # For other genes, check if the gene name is contained in the reference name
        elif gene_name.lower() in ref_gene.lower() or ref_gene.lower() in gene_name.lower():
            return color
    return "#808080"  # Default gray color for unmatched genes

##########################

def parse_hairpin_csv(file_path, regions, sequence_id):
    """
    Parse the hairpin CSV file to extract hairpin features for a specific sequence ID.
    """
    hairpins = []
    with open(file_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            # Check if the row corresponds to the requested sequence ID
            if not row['sequence_id'].lower().startswith(sequence_id.lower()):
                continue  # Skip rows that don't match the sequence ID

            # Skip rows with 'NA' in start_pos or end_pos
            if row['start_pos'] == 'NA' or row['end_pos'] == 'NA':
                continue

            region = row['region'].lower()
            start_pos = int(row['start_pos'])
            end_pos = int(row['end_pos'])
            
            # Adjust start and end positions based on the region
            if region == 'lsc':
                start_pos += regions['LSC'][0] - 1
                end_pos += regions['LSC'][0] - 1
            elif region == 'ssc':
                start_pos += regions['SSC'][0] - 1
                end_pos += regions['SSC'][0] - 1
            elif region == 'irb':
                start_pos += regions['IRb'][0] - 1
                end_pos += regions['IRb'][0] - 1
            
            hairpins.append({
                "id": row['ID'],
                "start": start_pos,
                "end": end_pos,
                "sequence": row['sequence'],
                "structure": row['structure'],
                "free_energy": float(row['free_energy']),
                "region": region,
                "genus": row['genus'],
                "species": row['species'],
                "sequence_id": row['sequence_id'],
                "major_clade": row['major_clade'],
                "subfamily": row['subfamily'],
                "tribe": row['tribe'],
                "subtribe": row['subtribe'],
                "photosynthetic_pathway": row['photosynthetic_pathway']
            })
    return hairpins

##########################

def reverse_genes_for_irb(ira_genes, ir_length):
    """
    Reverse the gene positions for IRb based on IRa genes.
    """
    irb_genes = []
    for gene in ira_genes:
        gene_name, (start, end) = gene
        new_start = ir_length - end + 1
        new_end = ir_length - start + 1
        irb_genes.append((gene_name, (new_start, new_end)))
    return irb_genes

##########################

def draw_genome_regions(ax, regions, region_colors):
    """
    Draw rectangles for genome regions on the plot.
    """
    ordered_regions = ['LSC', 'IRa', 'SSC', 'IRb']
    for region_name in ordered_regions:
        start, end = regions.get(region_name, (0, 0))
        color = region_colors.get(region_name, "#ffffff")  # Default to white if not found
        print(f"Drawing {region_name} from {start} to {end} with color {color}")  # Debugging output
        ax.add_patch(plt.Rectangle((start, -0.5), end - start, 1, color=color, alpha=0.8))

##########################

def plot_genome(regions, genes, tRNAs, hairpins, genome_length):
    """
    Create a linear plot of the chloroplast genome including hairpins.
    """
    fig, ax = plt.subplots(figsize=(15, 5))
    print(f"Genome length: {genome_length}. Regions: {regions}")
    ax.set_xlim(0, genome_length)
    ax.set_ylim(-1, 1)
    
    # Load gene colors from CSV
    gene_colors = load_gene_colors()
    
    # Define colors for genome regions
    region_colors = {
        "LSC": "#f5d893",
        "SSC": "#e8b26f",
        "IRa": "#b6834c",
        "IRb": "#b6834c"
    }
    
    # Draw genome regions
    draw_genome_regions(ax, regions, region_colors)
    
    count_pos = 0
    count_neg = 0
    # Plot genes
    for gene in genes:
        print(f"Gene: {gene}")
        if gene['start'] < gene['end'] and gene['start'] >= 0 and gene['end'] <= genome_length:
            y_pos = 0.6 if gene['strand'] == '+' else 0
            height = 0.4
            
            # Get color for this gene
            color = get_gene_color(gene['name'], gene_colors)

            ax.add_patch(Rectangle((gene["start"], y_pos), gene["end"] - gene["start"], height, 
                                 color=color, alpha=1))
            
            if (gene["end"] - gene["start"]) > 50:  # Only label genes above minimum size
                display_name = gene["id"] if gene["name"] == "Unknown" else gene["name"]
                text_y = 1.2 if gene['strand'] == '+' else -0.2
                ax.text((gene["start"] + gene["end"]) / 2, text_y, display_name, 
                       ha="center", va="center", fontsize=8, color=color,
                       rotation=30 if gene['strand'] == '+' else -20)

    # Plot tRNAs with simplified labeling
    for tRNA in tRNAs:
        # Check if a gene exists with the same start and end position
        print(f"tRNA: {tRNA}")
        gene_exists = any(gene['start'] == tRNA['start'] and gene['end'] == tRNA['end'] for gene in genes)
        
        if not gene_exists and tRNA['start'] < tRNA['end'] and tRNA['start'] >= 0 and tRNA['end'] <= genome_length:
            y_pos = 0.8 if tRNA['strand'] == '+' else -0.2
            height = 0.4
            
            ax.add_patch(Rectangle((tRNA["start"], y_pos), tRNA["end"] - tRNA["start"], height, 
                                 color="red", alpha=0.5))
            
            if (tRNA["end"] - tRNA["start"]) > 50:  # Only label tRNAs above minimum size
                display_name = tRNA["id"] if tRNA["name"] == "Unknown" else tRNA["name"]
                text_y = 1.8 if tRNA['strand'] == '+' else -0.8
                y_change = 0.16 if tRNA['strand'] == '+' else -0.16
            
                if tRNA['strand'] == '+':
                    count_pos += 1
                else:
                    count_neg += 1

                if count_pos == 2:
                    text_y = text_y + y_change

                if count_pos == 3:
                    text_y = text_y + (y_change * 2)
                    count_pos = 0
                # Draw a line from tRNA to label
                mid_x = (tRNA["start"] + tRNA["end"]) / 2
                ax.plot([mid_x, mid_x], [y_pos + height/2, text_y], linewidth=0.5)
                
                ax.text(mid_x, text_y, display_name, 
                       ha="center", va="center", fontsize=8,
                       rotation=30 if tRNA['strand'] == '+' else -30)

    # Plot hairpins
    for index, hairpin in enumerate(hairpins, start=1):
        #print(f"Hairpin: {hairpin}")
        y_center = 0.45
        height = 0.65
        ax.add_patch(Rectangle((hairpin["start"], y_center - height / 2), hairpin["end"] - hairpin["start"], height + 0.2, 
                             color="black", alpha=0.8, edgecolor="hotpink", linewidth=1))
        ax.add_patch(Rectangle((hairpin["start"], y_center - height / 2), hairpin["end"] - hairpin["start"], height, 
                             color="hotpink", alpha=1, edgecolor="hotpink", linewidth=1))
        mid_x = (hairpin["start"] + hairpin["end"]) / 2
        y_center = 0.2
        if index % 3 == 0:
            y_center = y_center - 0.3

        elif index % 2 == 0:
            y_center = y_center - 0.15

        ax.text(mid_x, y_center - height, hairpin["id"], 
                ha="center", va="center", fontsize=8, color="hotpink")

    # Add genome summary text
    summary_text = [
        f"LSC: {regions['LSC'][0]} - {regions['LSC'][1]} ({regions['LSC'][1] - regions['LSC'][0]} bp)",
        f"SSC: {regions['SSC'][0]} - {regions['SSC'][1]} ({regions['SSC'][1] - regions['SSC'][0]} bp)",
        f"IRa: {regions['IRa'][0]} - {regions['IRa'][1]} ({regions['IRa'][1] - regions['IRa'][0]} bp)",
        f"IRb: {regions['IRb'][0]} - {regions['IRb'][1]} ({regions['IRb'][1] - regions['IRb'][0]} bp)"
    ]
    num_hairpins = len(hairpins)  # Define num_hairpins as the length of the hairpins list
    summary_text.extend([
        f"Total Genes: {len(genes)}",
        f"Total tRNAs: {len(tRNAs)}",
        f"Total Hairpins: {num_hairpins}",
        f"Genome Length: {genome_length}"
    ])

    # Position the summary text in the upper left corner
    summary_y = 3.5
    summary_x = 0
    for line in summary_text:
        ax.text(summary_x, summary_y, line, 
                fontsize=8, va='top', ha='left')
        summary_y -= 0.15
        if summary_y <= 3.1:
            summary_y = 3.5
            summary_x += 25000

    gene_colors = {
        "ATP-synthase": "#4b726e",
        "Calvin-Cycle": "#4d4539",
        "Cytochrome-complex": "#8caba1",
        "Initiation-factor": "#79444a",
        "Large-subunit": "#ae5d40",
        "NADH-dehydrogenase": "#927441",
        "Other": "#847875",
        "Photosystem-I": "#77743b",
        "Photosystem-II": "#b3a555",
        "Small-subunit": "#c77b58",
        "Transcription-and-splicing": "#d1b187",
        "ribosomal-RNA": "#927441",
        "transfer-RNAs": "#ba9158"
    }
    
    # Print legend for each gene color
    for function, color in gene_colors.items():
        print(f"Function: {function}, Color: {color}")

    # Create legend handles for function-color pairs
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', label=function,
                                  markerfacecolor=color, markersize=10)
                       for function, color in gene_colors.items()]
    # Set plot limits and labels
    ax.set_xlim(0, genome_length)
    ax.set_ylim(-1, 3)  # Adjusted y-limits to accommodate labels
    ax.set_xlabel("Genome Position (bp)")
    ax.set_yticks([])
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1, 1), fontsize=8, ncol=5)

    # Set the title and subtitle using the first hairpin's metadata
    if hairpins:
        #print(f"Hairpins: {hairpins}")
        first_hairpin = hairpins[0]
        title = f"{first_hairpin['genus']} {first_hairpin['species']} ({first_hairpin['sequence_id']})"
        subtitle = f"{first_hairpin['major_clade']} | {first_hairpin['subfamily']} | {first_hairpin['tribe']} | {first_hairpin['subtribe']}. " \
                   f"Pathway = {first_hairpin['photosynthetic_pathway']}"
        ax.text(90000, 3.5, title, fontsize=14, va='top', ha='left', style='italic')
        ax.text(0.8, 1.05, subtitle, transform=ax.transAxes, fontsize=10, va='top', ha='center')

    plt.tight_layout()
    plt.show()


##########################

def filter_by_sequence_id(data, sequence_id):
    """
    Filter data by sequence ID, allowing partial and case-insensitive matches.
    """
    sequence_id_lower = sequence_id.lower()
    return [item for item in data if sequence_id_lower in item['sequence_id'].lower()]


##########################

def main(sequence_id, region_file=None, gff3_file=None, hairpin_file=None):
    """
    Main function to load, parse, and visualize the genome.
    """
    # Use default files if not provided
    if not region_file:
        region_file = "/Users/fdbramblepelt/Desktop/Lab/Hairpin_Chloroplast/Chloroplast_Vis/visualize_chloroplast_features/all.filtered_chloroplast_genome_region_lengths.gff3"
    if not gff3_file:
        gff3_file = "/Users/fdbramblepelt/Desktop/Lab/Hairpin_Chloroplast/gff3/file_bin/GFF3_genes_complete_genbank.csv"
    if not hairpin_file:
        hairpin_file = "/Users/fdbramblepelt/Desktop/Lab/Hairpin_Chloroplast/Chloroplast_Vis/visualize_chloroplast_features/trimmed_data.all_hairpins.csv"

    # Parse the regions file
    regions = parse_regions_file(region_file, sequence_id)

    # Calculate genome length
    genome_length = regions["LSC"][1] + regions["SSC"][1] + 2 * (regions["IR"][0][1] - regions["IR"][0][0])

    # Calculate region positions
    regions = calculate_region_positions(regions)

    # Debug output
    print(f"Regions after parsing: {regions}")

    # Parse the GFF3 file and filter by sequence ID
    features = parse_gff3_file(gff3_file, sequence_id)
    #features = filter_by_sequence_id(features, sequence_id)

    # Organize features into regions, genes, and tRNAs
    regions, genes, tRNAs = organize_features(features, regions)

    # Parse the hairpin file and filter by sequence ID
    hairpins = parse_hairpin_csv(hairpin_file, regions, sequence_id)
    hairpins = filter_by_sequence_id(hairpins, sequence_id)

    # Debug output
    print(f"Calculated genome length: {genome_length}")

    # Plot the genome
    plot_genome(regions, genes, tRNAs, hairpins, genome_length)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize chloroplast genome.")
    parser.add_argument("sequence_id", help="Sequence ID to search for in the input files.")
    parser.add_argument("--region_file", help="Path to the region file.", default=None)
    parser.add_argument("--gff3_file", help="Path to the GFF3 file.", default=None)
    parser.add_argument("--hairpin_file", help="Path to the hairpin CSV file.", default=None)

    args = parser.parse_args()
    main(args.sequence_id, args.region_file, args.gff3_file, args.hairpin_file)
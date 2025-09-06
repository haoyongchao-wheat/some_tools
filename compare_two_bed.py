#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import os
import sys
import argparse

def read_bed_file(file_path):
    """
    Read BED file and return DataFrame
    """
    columns = ['chromosome', 'start', 'end', 'gene_name', 'score', 'strand']
    df = pd.read_csv(file_path, sep='\t', header=None, names=columns)
    return df

def check_overlap(region1, region2):
    """
    Check if two genomic regions overlap
    region format: (chromosome, start, end)
    """
    chr1, start1, end1 = region1
    chr2, start2, end2 = region2
    
    # Must be on same chromosome
    if chr1 != chr2:
        return False
    
    # Check for overlap
    return not (end1 < start2 or end2 < start1)

def find_overlaps_and_unique(df1, df2):
    """
    Find overlapping regions and unique regions in each dataset
    """
    overlaps = []
    unique_method1 = []
    unique_method2 = []
    
    # Convert dataframes to lists of tuples for easier processing
    regions1 = [(row['chromosome'], row['start'], row['end'], row['gene_name'], row['strand']) 
                for _, row in df1.iterrows()]
    regions2 = [(row['chromosome'], row['start'], row['end'], row['gene_name'], row['strand']) 
                for _, row in df2.iterrows()]
    
    # Track which regions in df2 have been matched
    matched_in_df2 = set()
    
    # For each region in df1, check for overlaps in df2
    for i, region1 in enumerate(regions1):
        chr1, start1, end1, name1, strand1 = region1
        found_overlap = False
        
        for j, region2 in enumerate(regions2):
            chr2, start2, end2, name2, strand2 = region2
            
            if check_overlap((chr1, start1, end1), (chr2, start2, end2)):
                # Found overlap
                overlaps.append({
                    'chromosome': chr1,
                    'start': min(start1, start2),
                    'end': max(end1, end2),
                    'gene_name_method1': name1,
                    'gene_name_method2': name2,
                    'overlap_type': 'overlap',
                    'source': 'both_methods',
                    'strand_method1': strand1,
                    'strand_method2': strand2
                })
                matched_in_df2.add(j)
                found_overlap = True
                break
        
        if not found_overlap:
            # Region only in method1
            unique_method1.append({
                'chromosome': chr1,
                'start': start1,
                'end': end1,
                'gene_name': name1,
                'overlap_type': 'unique_method1',
                'source': 'NLR-Annotator_only',
                'strand': strand1
            })
    
    # Find regions only in method2
    for j, region2 in enumerate(regions2):
        if j not in matched_in_df2:
            chr2, start2, end2, name2, strand2 = region2
            unique_method2.append({
                'chromosome': chr2,
                'start': start2,
                'end': end2,
                'gene_name': name2,
                'overlap_type': 'unique_method2',
                'source': 'filter_gene_only',
                'strand': strand2
            })
    
    return overlaps, unique_method1, unique_method2

def create_output_bed(overlaps, unique_method1, unique_method2, output_file):
    """
    Create comprehensive BED file with all results
    """
    all_results = []
    
    # Add overlaps
    for overlap in overlaps:
        all_results.append([
            overlap['chromosome'],
            overlap['start'],
            overlap['end'],
            f"{overlap['gene_name_method1']}|{overlap['gene_name_method2']}",
            overlap['overlap_type'],
            overlap['source']
        ])
    
    # Add unique method1
    for unique in unique_method1:
        all_results.append([
            unique['chromosome'],
            unique['start'],
            unique['end'],
            unique['gene_name'],
            unique['overlap_type'],
            unique['source']
        ])
    
    # Add unique method2
    for unique in unique_method2:
        all_results.append([
            unique['chromosome'],
            unique['start'],
            unique['end'],
            unique['gene_name'],
            unique['overlap_type'],
            unique['source']
        ])
    
    # Sort by chromosome and start position
    all_results.sort(key=lambda x: (x[0], x[1]))
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write("#chromosome\tstart\tend\tgene_name\toverlap_status\tsource_info\n")
        for result in all_results:
            f.write('\t'.join(map(str, result)) + '\n')
    
    return all_results

def create_statistics_report(overlaps, unique_method1, unique_method2, df1, df2):
    """
    Create detailed statistics report
    """
    stats = {
        'total_method1': len(df1),
        'total_method2': len(df2),
        'overlapping_regions': len(overlaps),
        'unique_method1': len(unique_method1),
        'unique_method2': len(unique_method2)
    }
    
    # Calculate percentages
    stats['overlap_percentage_method1'] = (stats['overlapping_regions'] / stats['total_method1']) * 100
    stats['overlap_percentage_method2'] = (stats['overlapping_regions'] / stats['total_method2']) * 100
    
    # Chromosome distribution
    chr_stats_method1 = df1['chromosome'].value_counts().to_dict()
    chr_stats_method2 = df2['chromosome'].value_counts().to_dict()
    
    return stats, chr_stats_method1, chr_stats_method2

def create_visualization(stats, chr_stats_method1, chr_stats_method2, output_dir, output_prefix='nlr_comparison_analysis'):
    """
    Create visualization plots
    """
    plt.style.use('default')
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 12
    
    # Create figure with subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. Overall comparison pie chart
    labels = ['Overlapping', 'bed1 Only', 'bed2 Only']
    sizes = [stats['overlapping_regions'], stats['unique_bed1'], stats['unique_bed2']]
    colors = ['#2E8B57', '#FF6B6B', '#4ECDC4']
    
    ax1.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
    ax1.set_title('bed Comparison', fontweight='bold', fontsize=14)
    
    # 2. Method comparison bar chart
    methods = ['bed1', 'bed2']
    totals = [stats['total_method1'], stats['total_method2']]
    
    bars = ax2.bar(methods, totals, color=['#FF6B6B', '#4ECDC4'], alpha=0.7)
    ax2.set_title('Total Genes Detected by EachBed', fontweight='bold', fontsize=14)
    ax2.set_ylabel('Number of Genes')
    
    # Add value labels on bars
    for bar, total in zip(bars, totals):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 10,
                f'{total}', ha='center', va='bottom', fontweight='bold')
    
    # 3. Chromosome distribution for method 1
    chromosomes = sorted(chr_stats_method1.keys())
    counts1 = [chr_stats_method1[chr] for chr in chromosomes]
    
    ax3.bar(range(len(chromosomes)), counts1, color='#FF6B6B', alpha=0.7)
    ax3.set_title('Bed1: Genes per Chromosome', fontweight='bold', fontsize=14)
    ax3.set_xlabel('Chromosome')
    ax3.set_ylabel('Number of Genes')
    ax3.set_xticks(range(len(chromosomes)))
    ax3.set_xticklabels(chromosomes, rotation=45)
    
    # 4. Chromosome distribution for method 2
    counts2 = [chr_stats_method2.get(chr, 0) for chr in chromosomes]
    
    ax4.bar(range(len(chromosomes)), counts2, color='#4ECDC4', alpha=0.7)
    ax4.set_title('Bed2: Genes per Chromosome', fontweight='bold', fontsize=14)
    ax4.set_xlabel('Chromosome')
    ax4.set_ylabel('Number of Genes')
    ax4.set_xticks(range(len(chromosomes)))
    ax4.set_xticklabels(chromosomes, rotation=45)
    
    plt.tight_layout()
    
    # Save plots
    plt.savefig(os.path.join(output_dir, f'{output_prefix}.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, f'{output_prefix}.svg'), bbox_inches='tight')
    plt.show()

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Compare bed methods')
    parser.add_argument('method1_file', help='First BED file (e.g., bed1 results)')
    parser.add_argument('method2_file', help='Second BED file (e.g., bed2 results)')
    parser.add_argument('--output-prefix', '-o', default='test', 
                       help='Output file prefix (default: test)')
    
    args = parser.parse_args()
    
    # File paths
    method1_file = args.method1_file
    method2_file = args.method2_file
    output_bed_file = f'{args.output_prefix}_results.bed'
    output_stats_file = f'{args.output_prefix}_statistics.txt'
    
    print("Reading BED files...")
    # Read the two BED files
    df1 = read_bed_file(method1_file)  # bed1 results
    df2 = read_bed_file(method2_file)  # bed2 gene results
    
    print(f"Method 1 (bed1): {len(df1)} genes")
    print(f"Method 2 (bed2): {len(df2)} genes")
    
    print("\nAnalyzing overlaps and unique regions...")
    # Find overlaps and unique regions
    overlaps, unique_method1, unique_method2 = find_overlaps_and_unique(df1, df2)
    
    print(f"Overlapping regions: {len(overlaps)}")
    print(f"Unique to bed1: {len(unique_method1)}")
    print(f"Unique to bed2: {len(unique_method2)}")
    
    print("\nCreating output BED file...")
    # Create comprehensive output BED file
    all_results = create_output_bed(overlaps, unique_method1, unique_method2, output_bed_file)
    
    print("\nGenerating statistics...")
    # Generate statistics
    stats, chr_stats_method1, chr_stats_method2 = create_statistics_report(
        overlaps, unique_method1, unique_method2, df1, df2
    )
    
    # Write statistics to file
    with open(output_stats_file, 'w') as f:
        f.write("BED Comparison Statistics\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total genes detected by bed1: {stats['total_method1']}\n")
        f.write(f"Total genes detected by bed2: {stats['total_method2']}\n")
        f.write(f"Overlapping regions: {stats['overlapping_regions']}\n")
        f.write(f"Unique to bed1: {stats['unique_method1']}\n")
        f.write(f"Unique to bed2: {stats['unique_method2']}\n\n")
        f.write(f"Overlap percentage (vs bed1): {stats['overlap_percentage_method1']:.2f}%\n")
        f.write(f"Overlap percentage (vs bed2): {stats['overlap_percentage_method2']:.2f}%\n\n")
        
        f.write("Chromosome distribution - bed1:\n")
        for chr, count in sorted(chr_stats_method1.items()):
            f.write(f"  {chr}: {count} genes\n")
        
        f.write("\nChromosome distribution - bed2:\n")
        for chr, count in sorted(chr_stats_method2.items()):
            f.write(f"  {chr}: {count} genes\n")
    
    print("\nCreating visualizations...")
    # Create visualizations
    create_visualization(stats, chr_stats_method1, chr_stats_method2, '.', args.output_prefix)
    
    print(f"\nAnalysis complete!")
    print(f"Output files created:")
    print(f"  - {output_bed_file}: Comprehensive BED file with overlap annotations")
    print(f"  - {output_stats_file}: Detailed statistics report")
    print(f"  - {args.output_prefix}.png/svg: Visualization plots")

if __name__ == "__main__":
    main()

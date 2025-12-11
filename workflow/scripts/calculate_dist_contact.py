import cooltools
import cooler
import bioframe
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
import warnings
import pandas as pd
import sys
import argparse
from matplotlib.colors import LogNorm
import os

#================================
# function  
#=======================================

def load_data(clr_path, genome, n_proc):
    """
    Loads cooler file, genome data, and calculates expected contacts.
    """
    print(f"Loading cooler: {clr_path}")
    try:
        clr = cooler.Cooler(clr_path)
    except Exception as e:
        print(f"Error loading cooler file: {e}", file=sys.stderr)
        return None, None, None

    print(f"Fetching genome data for {genome}...")
    hg38_chromsizes = bioframe.fetch_chromsizes(genome)
    hg38_cens = bioframe.fetch_centromeres(genome)
    
    # Create chromosome arms view
    hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)
    
    # Filter arms to only those in the cooler file
    hg38_arms = hg38_arms[hg38_arms.chrom.isin(clr.chromnames)].reset_index(drop=True)

    print(f"Calculating expected-cis contacts with {n_proc} processes...")
    cvd = cooltools.expected_cis(
        clr=clr,
        view_df=hg38_arms,
        smooth=False,
        aggregate_smoothed=False,
        nproc=n_proc
    )

    cvd_filtered = cvd.loc[
        (cvd['dist_bp'] > 0) & (cvd['contact_frequency'] > 0)
    ].copy()
    
    print(f"Loaded and filtered data: {len(cvd_filtered)} valid contact pairs.")
    
    return clr, hg38_arms, cvd_filtered

def plot_per_arm_p(cvd_df, hg38_arms_df, out_file="P_s_per_arm.png"):
    """
    Plots the P(s) curve for each chromosome arm in a grid.
    """
    print(f"Plotting per-arm P(s) curves to {out_file}...")
    regions = hg38_arms_df['name'].unique()
    n_regions = len(regions)

    n_cols = int(np.ceil(np.sqrt(n_regions)))
    n_rows = int(np.ceil(n_regions / n_cols))

    f, axes = plt.subplots(
        n_rows, n_cols, 
        figsize=(n_cols * 3, n_rows * 3), # Dynamic figsize
        sharex=True, sharey=True
    )
    axes_flat = axes.flatten()

    for i, region in enumerate(regions):
        ax = axes_flat[i]
        region_data = cvd_df.loc[cvd_df['region1'] == region]
        
        ax.loglog(
            region_data['dist_bp'],
            region_data['contact_frequency'],
            '.', # Use dots for clarity
            ms=1, # Make markers small
            alpha=0.5
        )
        ax.set_title(region)
        ax.set_aspect(1.0)
        ax.grid(lw=0.5)

    axes_flat[0].set_ylim(bottom=10**-8)
    
    # Add shared labels to the whole figure
    f.supxlabel('separation, bp')
    f.supylabel('IC contact frequency')
    
    # Hide any unused axes
    for i in range(n_regions, len(axes_flat)):
        axes_flat[i].set_visible(False)
        
    plt.tight_layout()
    f.savefig(out_file, dpi=300)
    plt.close(f)
    print("...done.")


def plot_hexbin_p(cvd_df, out_file="P_s_hexbin.png"):
    """
    Plots an aggregated P(s) curve for all arms using a 2D hexbin.
    """
    print(f"Plotting aggregated hexbin P(s) curve to {out_file}...")
    f, ax = plt.subplots(1, 1, figsize=(7, 6))
    if cvd_df.empty:
            print("...WARNING: No valid contact pairs found. Writing empty plot.")
            ax.text(0.5, 0.5, "No Data", 
                    horizontalalignment='center', 
                    verticalalignment='center', 
                    transform=ax.transAxes,
                    fontsize=20,
                    color='gray')
    else:
        hb = ax.hexbin(
            cvd_df['dist_bp'],
            cvd_df['contact_frequency'],
            xscale='log',
            yscale='log',
            gridsize=100,
            norm=LogNorm(),
            cmap='inferno'
        )
        f.colorbar(hb, ax=ax, label='Count per bin')

    ax.set(
        xlabel='separation, bp',
        ylabel='IC contact frequency')
    ax.set_aspect(1.0)
    ax.grid(lw=0.5)
    ax.set_ylim(bottom=10**-8) 

    plt.tight_layout()
    f.savefig(out_file, dpi=300)
    plt.close(f)
    print("...done.")


def calculate_loglog_fits(cvd_df, hg38_arms_df, out_file="P_s_fits.csv"):
    """
    Performs a log-log linear fit for each region and all regions.
    Saves the results to a CSV file.
    """
    # (This function was already well-written, so I've left it as-is,
    #  but added a print statement and file saving)
    print("Calculating log-log linear fits...")
    # Define columns for the output table
    output_columns = ['region', 'slope', 'mse', 'n_points']
    
    # Check if the input DataFrame is empty
    if cvd_df.empty:
        print("...WARNING: No valid contact pairs found. Writing empty table.")
        
        # Create an empty DataFrame with the correct columns
        results_df = pd.DataFrame(columns=output_columns)
        
        # Add a placeholder "ALL_REGIONS" row
        all_region_stats = {
            'region': 'ALL_REGIONS', 'slope': np.nan, 
            'mse': np.nan, 'n_points': 0
        }
        all_region_df = pd.DataFrame([all_region_stats])
        results_df = pd.concat([results_df, all_region_df], ignore_index=True)

        # Save the empty table and exit the function
        results_df.to_csv(out_file, index=False)
        print(f"Results saved to {out_file}")
        return results_df    

    regions = hg38_arms_df['name'].unique()
    results_list = []
    
    warnings.filterwarnings('ignore', category=RuntimeWarning)

    for region in regions:
        region_data = cvd_df.loc[cvd_df['region1'] == region]
        
        # This filter is still good practice, even if cvd_df is pre-filtered,
        # as it makes the function "pure" (not relying on outside state)
        region_data = region_data.loc[
            (region_data['dist_bp'] > 0) & (region_data['contact_frequency'] > 0)
        ]
        
        if region_data.shape[0] < 2:
            results_list.append({
                'region': region, 'slope': np.nan, 
                'mse': np.nan, 'n_points': region_data.shape[0]
            })
            continue 

        region_data['log_dist_bp'] = np.log10(region_data['dist_bp'])
        region_data['log_contact_freq'] = np.log10(region_data['contact_frequency'])
        
        X = region_data[['log_dist_bp']]
        y = region_data['log_contact_freq']
        
        model = LinearRegression()
        model.fit(X, y)
        
        y_pred = model.predict(X)
        slope = model.coef_[0] 
        mse = mean_squared_error(y, y_pred)
        
        results_list.append({
            'region': region, 'slope': slope, 
            'mse': mse, 'n_points': len(y)
        })

    results_df = pd.DataFrame(results_list)
    
    # Calculate for ALL REGIONS
    # Note: cvd_df is already filtered, so this is fine
    if cvd_df.shape[0] >= 2:
        all_data = cvd_df.copy() # Use the pre-filtered dataframe
        all_data['log_dist_bp'] = np.log10(all_data['dist_bp'])
        all_data['log_contact_freq'] = np.log10(all_data['contact_frequency'])
        
        X_all = all_data[['log_dist_bp']]
        y_all = all_data['log_contact_freq']
        
        model_all = LinearRegression()
        model_all.fit(X_all, y_all)
        
        y_pred_all = model_all.predict(X_all)
        slope_all = model_all.coef_[0]
        mse_all = mean_squared_error(y_all, y_pred_all)
        
        all_region_stats = {
            'region': 'ALL_REGIONS', 'slope': slope_all, 
            'mse': mse_all, 'n_points': len(y_all)
        }
        all_region_df = pd.DataFrame([all_region_stats])
        results_df = pd.concat([results_df, all_region_df], ignore_index=True)

    print("...fits calculated.")
    
    ## CHANGE: Save the results
    results_df.to_csv(out_file, index=False)
    print(f"Results saved to {out_file}")
    
    return results_df

# -----------------------------------------------------------------
#  ARGUMENT PARSING
# -----------------------------------------------------------------

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Perform P(s) analysis on a Cooler file.")
    
    parser.add_argument(
        '-i', '--input',
        help="Input mcool file path (e.g., /path/to/file.mcool)",
        type=str,
        required=True
    )
    parser.add_argument(
        '-g', '--genome',
        help="Genome assembly name (e.g., hg38)",
        type=str,
        default='hg38'
    )
    parser.add_argument(
        '-r', '--resolution',
        help="Resolution (bin size) to use from the mcool file",
        type=int,
        default=1000
    )
    parser.add_argument(
        '-t', '--threads',
        help="Number of processes to use for calculation",
        type=int,
        default=40
    )
    parser.add_argument(
        '-p', '--prefix',
        help="Prefix for output files (e.g., 'my_analysis'). "
             "This will create 'my_analysis_per_arm.png', 'my_analysis_hexbin.png', etc.",
        type=str,
        required=True
    )
    parser.add_argument(
        '-o', '--output',
        help="Directory to output",
        type=str,
        required=True
    )
    
    return parser.parse_args()
# -----------------------------------------------------------------
#  MAIN EXECUTION
# -----------------------------------------------------------------

def main(args):
    """
    Main workflow for the P(s) analysis.
    """
    
    # 1. Construct paths from args
    cooler_uri = f"{args.input}::/resolutions/{args.resolution}"
    out_file_per_arm = f"{args.output}/{args.prefix}_per_arm.png"
    out_file_hexbin = f"{args.output}/{args.prefix}_hexbin_all_arms.png"
    out_file_fits = f"{args.output}/{args.prefix}_loglog_fits.csv"
    if not os.path.isdir(args.output):
        os.mkdir(args.output)
    # 2. Load and process data
    # Pass the cooler URI, genome, and thread count
    _, hg38_arms, cvd = load_data(cooler_uri, args.genome, args.threads)
    
    if cvd is None:
        print("Data loading failed. Exiting.", file=sys.stderr)
        sys.exit(1) # Exit with an error code

    # 3. Create plots, passing the constructed output paths
    plot_per_arm_p(cvd, hg38_arms, out_file=out_file_per_arm)
    plot_hexbin_p(cvd, out_file=out_file_hexbin)

    # 4. Perform analysis, passing the constructed output path
    res_table = calculate_loglog_fits(cvd, hg38_arms, out_file=out_file_fits)
    
    print("\n--- Analysis Complete ---")
    print("Slope/MSE Table (first 5 rows):")
    print(res_table.head())
    print("\nSlope/MSE for ALL_REGIONS:")
    print(res_table[res_table['region'] == 'ALL_REGIONS'])


if __name__ == "__main__":
    # 1. Parse arguments
    args = parse_args()
    # 2. Run the main workflow
    main(args)
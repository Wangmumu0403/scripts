import matplotlib.pyplot as plt
import numpy as np
import os

def plot_dp_test_results_scatter(output_files_config, save_path="dp_test.png"):
    """
    Plots DeePMD-kit test output files as scatter plots.
    Uses DFT (data) values on the X-axis and predicted values on the Y-axis.
    Skips the first header line in each file. Arranges plots in two rows.

    Args:
        output_files_config (list): A list of dictionaries, each defining a plot:
            e.g., {'filepath': 'test.e.out', 'title': 'Total Energy', 'ylabel': 'Energy (eV)', 'type': 'e'}
                  {'filepath': 'test.f.out', 'title': 'Forces', 'ylabel': 'Force (eV/Å)', 'type': 'f'}
                  {'filepath': 'test.v.out', 'title': 'Virial', 'ylabel': 'Virial (eV)', 'type': 'v'}
        save_path (str): The path to save the combined plot.
    """
    num_plots = len(output_files_config)
    # Arrange plots in 2 rows, calculate columns needed
    num_rows = 2
    num_cols = (num_plots + num_rows - 1) // num_rows # Ceiling division to get required columns

    # Adjust figure size for better readability, especially with multiple plots
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(6 * num_cols, 5 * num_rows))
    
    # Flatten axes array for easy iteration if it's 2D
    axes = axes.flatten()

    # Define component names for force (f) and virial (v) legends
    component_names_f = ["fx", "fy", "fz"]
    component_names_v = ["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"]

    for i, config in enumerate(output_files_config):
        filepath = config['filepath']
        title = config['title']
        ylabel_full = config['ylabel'] # Full Y-axis label, e.g., "Energy (eV)"
        plot_type = config['type']

        # Extract base name from full Y-axis label for X/Y axis titles
        ylabel_base = ylabel_full.split('(')[0].strip()

        ax = axes[i] # Get the current subplot axis

        if not os.path.exists(filepath):
            print(f"Warning: File not found: {filepath}. Skipping this plot.")
            ax.set_title(f"{title} (File not found)", color='red')
            # Clear axis if skipped, to prevent empty plot with default ticks
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        print(f"Processing file: {filepath} (Type: {plot_type})")
        try:
            # Load data, skipping the first row (header). unpack=True transposes the data.
            data_raw = np.loadtxt(filepath, skiprows=1, unpack=True) 

            ax.set_title(title, fontsize=14)
            ax.set_xlabel(f"DFT {ylabel_base}", fontsize=12)
            ax.set_ylabel(f"Predicted {ylabel_base}", fontsize=12)
            ax.grid(True, linestyle='--', alpha=0.7)

            min_val = np.inf
            max_val = -np.inf
            
            # Plot scatter based on file type
            if plot_type == 'e': # Energy files (test.e.out, test.e_peratom.out)
                # Data format: # dptest: data_e pred_e
                dft_vals = data_raw[0]
                pred_vals = data_raw[1]
                ax.scatter(dft_vals, pred_vals, s=5, alpha=0.7, label='Data points')
                
                min_val = min(dft_vals.min(), pred_vals.min())
                max_val = max(dft_vals.max(), pred_vals.max())
                
                # Add y=x reference line
                ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.6, label='y=x')
                ax.set_aspect('equal', adjustable='box') # Ensure square aspect ratio
                ax.legend(fontsize=10)

            elif plot_type == 'f': # Force file (test.f.out)
                # Data format: # dptest: data_fx data_fy data_fz pred_fx pred_fy pred_fz (6 columns)
                # data_raw[0]=data_fx, data_raw[1]=data_fy, data_raw[2]=data_fz
                # data_raw[3]=pred_fx, data_raw[4]=pred_fy, data_raw[5]=pred_fz
                
                # Plot scatter for each component
                for j in range(3): # For fx, fy, fz
                    dft_vals = data_raw[j]
                    pred_vals = data_raw[j+3]
                    ax.scatter(dft_vals, pred_vals, s=5, alpha=0.6, label=f'{component_names_f[j]}')
                    
                    # Update min/max for all components to set consistent axis limits
                    min_val = min(min_val, dft_vals.min(), pred_vals.min())
                    max_val = max(max_val, dft_vals.max(), pred_vals.max())
                
                # Add y=x reference line
                ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.6, label='y=x')
                ax.set_aspect('equal', adjustable='box')
                ax.legend(fontsize=10, loc='best', ncol=3) # Force component legend in 3 columns

            elif plot_type == 'v': # Virial files (test.v.out, test.v_peratom.out)
                # Data format: # dptest: data_vxx ... data_vzz pred_vxx ... pred_vzz (18 columns)
                # data_raw[0..8] are data_v components, data_raw[9..17] are pred_v components
                
                # Plot scatter for each tensor component
                for j in range(9): # For 9 tensor components
                    dft_vals = data_raw[j]
                    pred_vals = data_raw[j+9]
                    ax.scatter(dft_vals, pred_vals, s=5, alpha=0.5, label=f'{component_names_v[j]}')
                    
                    min_val = min(min_val, dft_vals.min(), pred_vals.min())
                    max_val = max(max_val, dft_vals.max(), pred_vals.max())
                
                # Add y=x reference line
                ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.6, label='y=x')
                ax.set_aspect('equal', adjustable='box')
                ax.legend(fontsize=8, loc='best', ncol=3) # Virial component legend in 3 columns

            # Set consistent axis limits and add a small buffer
            buffer = (max_val - min_val) * 0.05
            ax.set_xlim(min_val - buffer, max_val + buffer)
            ax.set_ylim(min_val - buffer, max_val + buffer)


        except Exception as e:
            print(f"Error processing {filepath}: {e}")
            print(f"Please check the content and format of {filepath}. It should contain correct numerical data. Error: {e}")
            # Clear axis if error occurs
            ax.set_title(f"{title} (Error)", color='red')
            ax.set_xticks([])
            ax.set_yticks([])
            continue
    
    # Remove any unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout() # Adjust layout to prevent overlapping titles/labels
    plt.savefig(save_path, dpi=300)
    print(f"\nPlot saved to {save_path}")
    plt.show()

# --- Configuration ---
if __name__ == "__main__":
    # Define your output files and their types. Adjust paths if they are not in the same directory.
    # Make sure these files exist before running the script.
    dp_output_configs = [
        {'filepath': "test.e.out", 'title': "Total Energy", 'ylabel': "Energy (eV)", 'type': 'e'},
        {'filepath': "test.e_peratom.out", 'title': "Energy per Atom", 'ylabel': "Energy/atom (eV)", 'type': 'e'},
        {'filepath': "test.f.out", 'title': "Force Components", 'ylabel': "Force (eV/Å)", 'type': 'f'},
        {'filepath': "test.v.out", 'title': "Virial Tensor", 'ylabel': "Virial (eV)", 'type': 'v'},
        {'filepath': "test.v_peratom.out", 'title': "Virial Tensor per Atom", 'ylabel': "Virial/atom (eV)", 'type': 'v'}
    ]

    # Output path for the combined plot
    plot_output_filename = "dp_test.png"

    # --- Run the plotting function ---
    plot_dp_test_results_scatter(dp_output_configs, plot_output_filename)

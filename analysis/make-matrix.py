import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.polynomial import polyfit
import re
import matplotlib as mpl

from utility import get_all_stats, fmt


def make_publication_plot2(dfs,
                           figsize=(20, 18),
                           replicates=100_000,
                           confidence=95, 
                           min_value=-12,
                           max_value=0,
                           pad=1):

    mpl.rcParams['font.size'] = 28
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'  # Use Computer Modern fonts
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath,amsfonts}'

    def get_numerical_values(text):
        pattern = r"(-?\d+\.\d+)\s+\[\d+%:\s+(-?\d+\.\d+),\s+(-?\d+\.\d+)\]"
        match = re.search(pattern, text)
        if match:
            value = float(match.group(1))  # First number
            min_value = float(match.group(2))  # Min of the interval
            max_value = float(match.group(3))  # Max of the interval
            return value, min_value, max_value

    FF_NAME = {
        'espaloma-0.3.1': r'\textbf{Espaloma}',
        'espaloma-0.3.1*': r'\textbf{Espaloma*}',
        'gaff-2.11': r'\textbf{GAFF}',
        'openff-2.0.0': r'\textbf{OpenFF}',
        'Chen2023': r'\textbf{Chen et al. 2023}'
    }

    # Determine the range for the plot
    simulations = [column for column in dfs[0][1].columns if column.startswith('simulation')]

    # Create the plot
    fig, axes = plt.subplots(nrows=len(dfs), ncols=len(simulations), figsize=figsize, sharex=True, sharey=True)

    for i, item in enumerate(dfs):
        system_name, df = item
        system_name = SYSTEM_NAME[system_name]
        # Get all stats of the simulations
        all_stats = get_all_stats(df, replicates=replicates, confidence=confidence)

        for j, simulation in enumerate(simulations):
            filtered_df = df[df[simulation].notna()]
            if len(filtered_df):
                axes[i, j].grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

                pearson = fmt(get_numerical_values(all_stats.loc[simulation, 'pearson']))
                spearman = fmt(get_numerical_values(all_stats.loc[simulation, 'spearman']))
                kendall = fmt(get_numerical_values(all_stats.loc[simulation, 'kendall']))
                rmse = fmt(get_numerical_values(all_stats.loc[simulation, 'rmse']))
                mse = fmt(get_numerical_values(all_stats.loc[simulation, 'mse']))
                mue = fmt(get_numerical_values(all_stats.loc[simulation, 'mue']))

                sys_info_text = (
                    r"~\\"
                    fr"$\mathrm{{N}} = {{{len(df)}}}$"
                )

                stats_text = (
                    r"\begin{equation*}"
                    r"\begin{aligned}[t]"
                    fr"  &\rho = {pearson[0]}_{{{pearson[1]}}}^{{{pearson[2]}}} \\"
                    fr"  &r_S = {spearman[0]}_{{{spearman[1]}}}^{{{spearman[2]}}} \\"
                    fr"  &\tau = {kendall[0]}_{{{kendall[1]}}}^{{{kendall[2]}}} \\"
                    fr"  &\mathrm{{MSE}} = {mse[0]}_{{{mse[1]}}}^{{{mse[2]}}} \\"
                    fr"  &\mathrm{{MUE}} = {mue[0]}_{{{mue[1]}}}^{{{mue[2]}}} \\"
                    fr"  &\mathrm{{RMSE}} = {rmse[0]}_{{{rmse[1]}}}^{{{rmse[2]}}}"
                    r"\end{aligned}"
                    r"\end{equation*}"
                )

                axes[i, j].text(
                    0.02, 0.98, stats_text, transform=axes[i, j].transAxes, fontsize=19,
                    verticalalignment='top', horizontalalignment='left',
                    # bbox=dict(facecolor='white', alpha=0.9, edgecolor='gray')
                    )

                axes[i, j].text(
                    0.98, 0.02, sys_info_text, transform=axes[i, j].transAxes, fontsize=19,
                    verticalalignment='bottom', horizontalalignment='right',
                    )

                simulation_name = simulation.split('_')[-1]
                # Plot the y = x line
                axes[i, j].plot([min_value, max_value], [min_value, max_value], '-', label='Y = X', color='black')
                # Fill the region between the error curves
                error = 2
                axes[i, j].fill_between([min_value, max_value], [min_value - error, max_value - error], [min_value + error, max_value + error],
                                        alpha=0.3, label=r'$\pm$ 2 kcal/mol', color='gray')
                error = 1
                axes[i, j].fill_between([min_value, max_value], [min_value - error, max_value - error], [min_value + error, max_value + error],
                                        alpha=0.3, label=r'$\pm$ 1 kcal/mol', color='gray')

                b, m = polyfit(filtered_df.exp_dG, filtered_df[simulation], 1)
                inter_range = np.array([min_value, max_value])
                axes[i, j].plot(inter_range, b + m * inter_range, '--', label='Correlation line', color='black')

                # Plot the data

                # Scatter plot with color scale
                scatter = axes[i, j].scatter(
                    filtered_df.exp_dG, filtered_df[simulation],
                    c=(filtered_df.exp_dG - filtered_df[simulation]).abs(),
                    cmap=plt.cm.coolwarm,
                    norm=plt.Normalize(0, 5),
                    s=80, edgecolor='black', zorder=2)  # Higher zorder for dots

                # Add error bars behind the scatter points
                axes[i, j].errorbar(
                    filtered_df.exp_dG, filtered_df[simulation],
                    fmt='none', ecolor='gray',
                    xerr=filtered_df.exp_dG_error, yerr=filtered_df[f"sem_{simulation_name}"], zorder=1)

                if i == 0:
                    axes[i, j].set_title(FF_NAME[simulation_name], pad=20)

                # Set custom x- and y-axis labels with the standard-state notation
                if j == 0:
                    axes[i, j].set_ylabel(
                        fr"$\mathrm{{\textbf{{{system_name}}}}}$"+"\n"+r"$\Delta G_{{\mathrm{{calc}}}} \, \mathrm{{[kcal/mol]}}$", 
                        fontsize=28,
                        ha='center',  # horizontal alignment
                    )

                if i == len(dfs) - 1:
                    axes[i, j].set_xlabel(r"$\Delta G_{\mathrm{exp}} \mathrm{[kcal/mol]}$", fontsize=28)

                # Set equal scaling for x and y axes
                axes[i, j].set_aspect('equal', adjustable='box')

                axes[i, j].set_xlim([min_value, max_value])
                axes[i, j].set_ylim([min_value, max_value])

                # Compute the start and end values after applying padding
                start = min_value + pad
                end = max_value - pad/2

                # Calculate the step size for four equidistant points
                step = int((end - start) / 3)  # Divide the range into 3 steps for 4 points

                # Generate the ticks as integers
                start = round(start)
                ticks = [
                    start,
                    start + step,
                    start + 2 * step,
                    start + 3 * step,
                ]

                axes[i, j].set_xticks(ticks)
                axes[i, j].set_yticks(ticks)
            else:
                axes[i, j].axis('off')  # Turn off the axis for empty columns

                axes[i-1, j].set_xlabel(r"$\Delta G_{\mathrm{exp}} \mathrm{[kcal/mol]}$", fontsize=28)
                # Define the tick values based on the common limits
                # Compute the start and end values after applying padding
                start = min_value + pad
                end = max_value - pad/2
                # Calculate the step size for four equidistant points
                step = int((end - start) / 3)  # Divide the range into 3 steps for 4 points
                # Generate the ticks as integers
                start = round(start)
                ticks = [
                    start,
                    start + step,
                    start + 2 * step,
                    start + 3 * step,
                ]

                axes[i-1, j].set_xticks(ticks)
                axes[i-1, j].tick_params(axis='x', labelbottom=True)

    fig.subplots_adjust(wspace=0.01, hspace=0.1)  # Adjust the value as needed

    # Add colorbar
    cbar = fig.colorbar(scatter, ax=axes, location="right", shrink=0.8, pad=0.01, aspect=40)  # Attach the colorbar to the current axis
    cbar.set_label(r'$\lvert \Delta G_{\mathrm{calc}} -  \Delta G_{\mathrm{exp}} \rvert$ [kcal/mol]', fontsize=28)
    return axes, fig


if __name__ == "__main__":

    replicates = 10_000
    confidence = 68

    # Define the custom orders
    SYSTEM_ORDER = ["p38", "A2A", "ptp1b", "tyk2", "thrombin", "mcl1", "CyclophilinD", "SAMPL6-OA"]
    FORCE_FIELD_ORDER = ['espaloma-0.3.1', 'gaff-2.11', 'openff-2.0.0']

    SYSTEM_NAME = {
        "p38": "P38", 
        "A2A": "A2A",
        "ptp1b": "PTP1B",
        "tyk2": "TYK2",
        "thrombin": "Thrombin",
        "mcl1": "MCL1",
        "CyclophilinD": "CyclophilinD",
        "SAMPL6-OA": "SAMPL6-OA"
    }

    for CALC_TYPE in ["mbar", "dh-gb"]:
        BindFlowData = pd.read_csv("../data/simulation/bindflow/gather/BindFlow.csv", index_col=0)

        if CALC_TYPE == "dh-gb":
            pad = 10
        else:
            pad = 1

        columns = [
            "system",
            "ligand",
            "replica",
            "sample",
            "exp_dG",
            "exp_dG_error",
            f"simulation_{CALC_TYPE}_espaloma-0.3.1",
            f"simulation_{CALC_TYPE}_gaff-2.11",
            f"simulation_{CALC_TYPE}_openff-2.0.0",
        ]
        BindFlowData = BindFlowData[columns]

        BindFlowData.rename(
            columns={
                "system": "source",
                f"simulation_{CALC_TYPE}_espaloma-0.3.1": "simulation_espaloma-0.3.1",
                f"simulation_{CALC_TYPE}_gaff-2.11": "simulation_gaff-2.11",
                f"simulation_{CALC_TYPE}_openff-2.0.0": "simulation_openff-2.0.0",
            },
            inplace=True
        )

        mean = BindFlowData.groupby(["source", "ligand"]).mean().reset_index().drop(columns=["replica", "sample"])
        sem = BindFlowData.groupby(["source", "ligand"]).sem().reset_index().drop(columns=["replica", "sample"])


        # Filter DataFrame
        mask = mean["simulation_espaloma-0.3.1"].notna()
        mean = mean[mask]
        sem = sem[mask]

        sem.rename(
            columns={
                "simulation_espaloma-0.3.1": "sem_espaloma-0.3.1",
                "simulation_gaff-2.11": "sem_gaff-2.11",
                "simulation_openff-2.0.0": "sem_openff-2.0.0",
            },
            inplace=True
        )
        sem.drop(columns=["exp_dG", "exp_dG_error"], inplace=True)
        df_merge = pd.merge(mean, sem, on=["source", "ligand"])


        dfs1 = []
        for system in SYSTEM_ORDER[:4]:
            dfs1.append((system, df_merge[df_merge["source"] == system]))

        dfs2 = []
        for system in SYSTEM_ORDER[4:]:
            dfs2.append((system, df_merge[df_merge["source"] == system]))

        min_value = np.min(df_merge[["exp_dG"] + [f"simulation_{ff}" for ff in FORCE_FIELD_ORDER]]) - 1
        max_value = np.max(df_merge[["exp_dG"] + [f"simulation_{ff}" for ff in FORCE_FIELD_ORDER]]) + 1

        print(CALC_TYPE)
        ax, fig = make_publication_plot2(dfs1, figsize=(22, 24.5), replicates=replicates, confidence=confidence,
                                         min_value=min_value, max_value=max_value, pad=pad)
        fig.savefig(f"summary/matrix/{CALC_TYPE}-matrix1.pdf",
                    bbox_inches="tight",  # trims whitespace
                    pad_inches=0.0,       # removes extra padding
                    transparent=False)

        ax, fig = make_publication_plot2(dfs2, figsize=(22, 24.5), replicates=replicates, confidence=confidence,
                                         min_value=min_value, max_value=max_value, pad=pad)
        fig.savefig(f"summary/matrix/{CALC_TYPE}-matrix2.pdf",
                    bbox_inches="tight",  # trims whitespace
                    pad_inches=0.0,       # removes extra padding
                    transparent=False)

        # Warning are coming from attending to compute SAMPL6-OA with openff and gaff (data does not exists)

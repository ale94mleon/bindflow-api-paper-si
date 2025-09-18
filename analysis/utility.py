import re
from concurrent.futures import ProcessPoolExecutor
from typing import Callable, Iterable, List, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython.display import display
from numpy.polynomial.polynomial import polyfit
from plbenchmark import edges, ligands, targets


def get_numerical_values(text):
    pattern = r"(-?\d+\.\d+)\s+\[\d+%:\s+(-?\d+\.\d+),\s+(-?\d+\.\d+)\]"
    match = re.search(pattern, text)
    if match:
        value = float(match.group(1))  # First number
        min_value = float(match.group(2))  # Min of the interval
        max_value = float(match.group(3))  # Max of the interval
        return [value, min_value, max_value]
    else:
        return [None, None, None]

def fmt(values):
    """Format main, lower, upper values with 2 decimals."""
    return tuple(f"{v:.2f}" for v in values)

def make_publication_plot(df, figsize=(20, 16), replicates=100_000, confidence=95, make_scale=False, scale_on_grpah=False):
    import re

    import matplotlib as mpl

    mpl.rcParams['font.size'] = 28
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'  # Use Computer Modern fonts
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath,amsfonts}'


    # Determine the range for the plot
    simulations = [column for column in df.columns if column.startswith('simulation')]
    min_value = int(min(np.min(df.exp_dG), min([np.min(df[simulation]) for simulation in simulations])) - 1)
    max_value = int(max(np.max(df.exp_dG), max([np.max(df[simulation]) for simulation in simulations])) + 1)
    # Create the plot

    fig, axes = plt.subplots(nrows=1, ncols=len(simulations), figsize=figsize)

    if len(simulations) > 1:
        axes = axes.flatten()
    else:
        axes = [axes]

    # Get all stats of the simulations
    all_stats = get_all_stats(df, replicates=replicates, confidence=confidence)

    i = 0
    for simulation, ax in zip(simulations, axes):
        filtered_df = df[df[simulation].notna()]

        ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

        stats_text = (
            r"\textbf{Statistics:} \\[0.5em] "
            fr"$\mathrm{{Pearson}}={get_numerical_values(all_stats.loc[simulation, 'pearson'])[0]}_{{{get_numerical_values(all_stats.loc[simulation, 'pearson'])[1]}}}^{{{get_numerical_values(all_stats.loc[simulation, 'pearson'])[2]}}}$ \\[0.3em] "
            fr"$\mathrm{{Kendall}}={get_numerical_values(all_stats.loc[simulation, 'kendall'])[0]}_{{{get_numerical_values(all_stats.loc[simulation, 'kendall'])[1]}}}^{{{get_numerical_values(all_stats.loc[simulation, 'kendall'])[2]}}}$ \\[0.3em] "
            fr"$\mathrm{{Spearman}}={get_numerical_values(all_stats.loc[simulation, 'spearman'])[0]}_{{{get_numerical_values(all_stats.loc[simulation, 'spearman'])[1]}}}^{{{get_numerical_values(all_stats.loc[simulation, 'spearman'])[2]}}}$ \\[0.3em] "
            # r"Predictive Index: $0.44_{0.06}^{0.71}$ \\ "
            fr"$\mathrm{{RMSE}}={get_numerical_values(all_stats.loc[simulation, 'rmse'])[0]}_{{{get_numerical_values(all_stats.loc[simulation, 'rmse'])[1]}}}^{{{get_numerical_values(all_stats.loc[simulation, 'rmse'])[2]}}} \mathrm{{~kcal/mol}}$ \\[0.3em] "
            fr"$\mathrm{{MSE}}={get_numerical_values(all_stats.loc[simulation, 'mse'])[0]}_{{{get_numerical_values(all_stats.loc[simulation, 'mse'])[1]}}}^{{{get_numerical_values(all_stats.loc[simulation, 'mse'])[2]}}} \mathrm{{~kcal/mol}}$\\[0.3em] "
            fr"$\mathrm{{MUE}}={get_numerical_values(all_stats.loc[simulation, 'mue'])[0]}_{{{get_numerical_values(all_stats.loc[simulation, 'mue'])[1]}}}^{{{get_numerical_values(all_stats.loc[simulation, 'mue'])[2]}}} \mathrm{{~kcal/mol}}$"
            )

        ax.text(
            0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=16,
            verticalalignment='top', horizontalalignment='left',
            bbox=dict(facecolor='white', alpha=0.9, edgecolor='gray')
            )

        simulation_name = simulation.split('_')[-1]
        # Plot the y = x line
        ax.plot([min_value, max_value], [min_value, max_value], '-', label='Y = X', color='black')
        # Fill the region between the error curves
        error = 2
        ax.fill_between([min_value, max_value], [min_value - error, max_value - error], [min_value + error, max_value + error],
                        alpha=0.3, label=r'$\pm$ 2 kcal/mol', color='gray')
        error = 1
        ax.fill_between([min_value, max_value], [min_value - error, max_value - error], [min_value + error, max_value + error],
                        alpha=0.3, label=r'$\pm$ 1 kcal/mol', color='gray')

        b, m = polyfit(filtered_df.exp_dG, filtered_df[simulation], 1)
        inter_range = np.array([min_value, max_value])
        ax.plot(inter_range, b + m * inter_range, '--', label='Correlation line', color='black')

        # Plot the data

        # Scatter plot with color scale
        scatter = ax.scatter(
            filtered_df.exp_dG, filtered_df[simulation],
            c=(filtered_df.exp_dG - filtered_df[simulation]).abs(),
            cmap=plt.cm.coolwarm,
            norm=plt.Normalize(0, 5),
            s=100, edgecolor='black', zorder=2)  # Higher zorder for dots

        # Add error bars behind the scatter points
        ax.errorbar(
            filtered_df.exp_dG, filtered_df[simulation],
            fmt='none', ecolor='gray',
            xerr=filtered_df.exp_dG_error, yerr=filtered_df[f"sem_{simulation_name}"], zorder=1)

        ax.set_title(simulation_name)

        # Set custom x- and y-axis labels with the standard-state notation
        if i == 0:
            ax.set_ylabel(r"$\Delta G_{\mathrm{calc}} \mathrm{[kcal/mol]}$", fontsize=28)
        ax.set_xlabel(r"$\Delta G_{\mathrm{exp}} \mathrm{[kcal/mol]}$", fontsize=28)


        # Set equal scaling for x and y axes
        ax.set_aspect('equal', adjustable='box')

        ax.set_xlim([min_value, max_value])
        ax.set_ylim([min_value, max_value])

        # # Define the tick values based on the common limits
        ticks = [min_value, int(min_value + (max_value - min_value)/4), int(min_value + (max_value - min_value)/2), int(min_value + 3*(max_value - min_value)/4), max_value]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)

        # ax.legend(fontsize=16)
        if i > 0:
            ax.set_yticklabels([])
        i += 1
    fig.subplots_adjust(wspace=0.1)  # Adjust the value as needed
    # fig.tight_layout()

    if make_scale:
        if scale_on_grpah:
            # Add colorbar
            cbar = fig.colorbar(scatter, ax=axes, location="right", shrink=0.5, pad=0.05)  # Attach the colorbar to the current axis
            cbar.set_label(r'$\lvert \Delta G_{\mathrm{exp}} - \Delta G_{\mathrm{calc}} \rvert$ [kcal/mol]', fontsize=28)
            return ax, fig
        else:
            # Create a new figure for the colorbar (horizontal style)
            cbar_fig = plt.figure(figsize=(2*figsize[0]/3, figsize[1]))  # Adjust the size (width, height) as needed
            # Create an axis for the colorbar at the bottom of the figure
            cbar_ax = cbar_fig.add_axes([0.1, 0.4, 0.8, 0.05])  # Adjust position (x, y, width, height)
            # Create a horizontal colorbar using the same colormap as your scatter plot
            cbar = plt.colorbar(scatter, cax=cbar_ax, orientation='horizontal')
            # Set the label for the colorbar
            cbar.set_label(r'$\lvert \Delta G_{\mathrm{exp}} - \Delta G_{\mathrm{calc}} \rvert$ [kcal/mol]', fontsize=28)

            # Show the colorbar figure
            # cbar_fig.tight_layout()
            return ax, fig, cbar_ax, cbar_fig
    else:
        return ax, fig




def make_publication_plot2(dfs, name_list, figsize=(20, 18), replicates=100_000, confidence=95):

    import matplotlib as mpl

    mpl.rcParams['font.size'] = 28
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'  # Use Computer Modern fonts
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath,amsfonts}'

    FF_NAME = {
        'espaloma-0.3.1': r'\textbf{Espaloma}',
        'espaloma-0.3.1*': r'\textbf{Espaloma*}',
        'gaff-2.11': r'\textbf{GAFF}',
        'openff-2.0.0': r'\textbf{OpenFF}',
        'Chen2023': r'\textbf{Chen et al. 2023}'
    }

    # Determine the range for the plot
    simulations = [column for column in dfs[0].columns if column.startswith('simulation')]
    min_value = np.inf
    max_value = -1*np.inf
    for df in dfs:
        min_value = int(min(np.min(df[simulations+["exp_dG"]]), min_value))
        max_value = int(max(np.max(df[simulations+["exp_dG"]]), min_value))
    min_value -= 1
    max_value += 1

    # Create the plot
    fig, axes = plt.subplots(nrows=len(dfs), ncols=len(simulations), figsize=figsize, sharex=True, sharey=True)

    for i, df in enumerate(dfs):
        # Get all stats of the simulations
        all_stats = get_all_stats(df, replicates=replicates, confidence=confidence)

        for j, simulation in enumerate(simulations):
            filtered_df = df[df[simulation].notna()]

            axes[i, j].grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

            sys_info_text = (
                r"~\\"
                fr"$\mathrm{{N}} = {{{len(df)}}}$"
            )
            if i == len(name_list) - 1:
                pearson = fmt(get_numerical_values(all_stats.loc[simulation, 'pearson']))
                spearman = fmt(get_numerical_values(all_stats.loc[simulation, 'spearman']))
                kendall = fmt(get_numerical_values(all_stats.loc[simulation, 'kendall']))
                rmse = fmt(get_numerical_values(all_stats.loc[simulation, 'rmse']))
                mue = fmt(get_numerical_values(all_stats.loc[simulation, 'mue']))

                stats_text = (
                    r"\begin{equation*}"
                    r"\begin{aligned}[t]"
                    fr"  &\rho = {pearson[0]}_{{{pearson[1]}}}^{{{pearson[2]}}} \\"
                    fr"  &r_S = {spearman[0]}_{{{spearman[1]}}}^{{{spearman[2]}}} \\"
                    fr"  &\tau = {kendall[0]}_{{{kendall[1]}}}^{{{kendall[2]}}} \\"
                    fr"  &\mathrm{{ocRMSE}} = {rmse[0]}_{{{rmse[1]}}}^{{{rmse[2]}}} \\"
                    fr"  &\mathrm{{ocMUE}} = {mue[0]}_{{{mue[1]}}}^{{{mue[2]}}}"
                    r"\end{aligned}"
                    r"\end{equation*}"
                )

            else:
                pearson = fmt(get_numerical_values(all_stats.loc[simulation, 'pearson']))
                spearman = fmt(get_numerical_values(all_stats.loc[simulation, 'spearman']))
                kendall = fmt(get_numerical_values(all_stats.loc[simulation, 'kendall']))
                rmse = fmt(get_numerical_values(all_stats.loc[simulation, 'rmse']))
                mse = fmt(get_numerical_values(all_stats.loc[simulation, 'mse']))
                mue = fmt(get_numerical_values(all_stats.loc[simulation, 'mue']))

                stats_text = (
                    r"\begin{equation*}"
                    r"\begin{aligned}[t]"
                    fr"  &\rho = {pearson[0]}_{{{pearson[1]}}}^{{{pearson[2]}}} \\"
                    fr"  &r_S = {spearman[0]}_{{{spearman[1]}}}^{{{spearman[2]}}} \\"
                    fr"  &\tau = {kendall[0]}_{{{kendall[1]}}}^{{{kendall[2]}}} \\"
                    fr"  &\mathrm{{RMSE}} = {rmse[0]}_{{{rmse[1]}}}^{{{rmse[2]}}} \\"
                    fr"  &\mathrm{{MSE}} = {mse[0]}_{{{mse[1]}}}^{{{mse[2]}}} \\"
                    fr"  &\mathrm{{MUE}} = {mue[0]}_{{{mue[1]}}}^{{{mue[2]}}}"
                    r"\end{aligned}"
                    r"\end{equation*}"
                )


            axes[i, j].text(
                0.98, 0.02, sys_info_text, transform=axes[i, j].transAxes, fontsize=22,
                verticalalignment='bottom', horizontalalignment='right',
                )
            axes[i, j].text(
                0.02, 0.98, stats_text, transform=axes[i, j].transAxes, fontsize=20,
                verticalalignment='top', horizontalalignment='left',
                # bbox=dict(facecolor='white', alpha=0.9, edgecolor='gray')
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
                c=(filtered_df[simulation] - filtered_df.exp_dG).abs(),
                cmap=plt.cm.coolwarm,
                norm=plt.Normalize(0, 5),
                s=100, edgecolor='black', zorder=2)  # Higher zorder for dots

            # Add error bars behind the scatter points
            axes[i, j].errorbar(
                filtered_df.exp_dG, filtered_df[simulation],
                fmt='none', ecolor='gray',
                xerr=filtered_df.exp_dG_error, yerr=filtered_df[f"sem_{simulation_name}"], zorder=1)

            if i == 0:
                axes[i, j].set_title(FF_NAME[simulation_name], pad=20)

            # Set custom x- and y-axis labels with the standard-state notation
            if j == 0:
                if i == len(dfs) - 1:
                    axes[i, j].set_ylabel(
                        fr"$\mathrm{{\textbf{{{name_list[i]}}}}}$"+"\n"+r"$\Delta G^{{\mathrm{{oc}}}}_{{\mathrm{{calc}}}} \, \mathrm{{[kcal/mol]}}$",
                        fontsize=28,
                        ha='center',  # horizontal alignment
                    )
                else:
                    axes[i, j].set_ylabel(
                        fr"$\mathrm{{\textbf{{{name_list[i]}}}}}$"+"\n"+r"$\Delta G_{{\mathrm{{calc}}}} \, \mathrm{{[kcal/mol]}}$",
                        fontsize=28,
                        ha='center',  # horizontal alignment
                    )

            if i == len(dfs) - 1:
                axes[i, j].set_xlabel(r"$\Delta G_{\mathrm{exp}} \mathrm{[kcal/mol]}$", fontsize=28)

            # Set equal scaling for x and y axes
            axes[i, j].set_aspect('equal', adjustable='box')

            axes[i, j].set_xlim([min_value, max_value])
            axes[i, j].set_ylim([min_value, max_value])

            # Define the tick values based on the common limits
            pad = 1
            # Compute the start and end values after applying padding
            start = (min_value + pad)
            end = max_value - pad/2
            # Calculate the step size for four equidistant points
            step = int((end - start) / 3)  # Divide the range into 3 steps for 4 points
            start = round(start)
            # Generate the ticks as integers
            ticks = [
                start,
                start + step,
                start + 2 * step,
                start + 3 * step,
            ]

            axes[i, j].set_xticks(ticks)
            axes[i, j].set_yticks(ticks)

        fig.subplots_adjust(wspace=0.01, hspace=0.1)  # Adjust the value as needed

        # Add colorbar
        cbar = fig.colorbar(scatter, ax=axes, location="right", shrink=0.8, pad=0.01)  # Attach the colorbar to the current axis
        cbar.set_label(r'$\lvert \Delta G_{\mathrm{calc}} ~(\mathrm{or}~\Delta G^{{\mathrm{{oc}}}}_{{\mathrm{{calc}}}}) - \Delta G_{\mathrm{exp}} \rvert$ [kcal/mol]', fontsize=28)
    return axes, fig


def get_all_stats(df, replicates=1000, confidence=95):
    all_stats = pd.merge(
        get_corr(df, replicates=replicates, confidence=confidence),
        get_rmse(df, replicates=replicates, confidence=confidence),
        left_index=True,
        right_index=True)
    return all_stats


# def get_corr(df, replicates=1000, confidence=95):
#     tail = (1 - confidence / 100) / 2
#     values = dict()

#     for method in ['pearson', 'kendall', 'spearman', 'predictive_index']:
#         values[method] = dict()
#         for column in df.columns:
#             if column.startswith('simulation_'):
#                 # Filter out rows with None or NaN in 'exp_dG' or the current simulation column
#                 filtered_df = df[['exp_dG', column]].dropna()

#                 # Estimate confidence intervals
#                 bootstrapped_values = []
#                 for indexes in get_bootstrap_indexes(filtered_df, replicates=replicates):
#                     bt_df = filtered_df[['exp_dG', column]].iloc[indexes, :]
#                     if method == 'predictive_index':
#                         bootstrapped_values.append(bt_df.corr(method=predictive_index).loc['exp_dG', column])
#                     else:
#                         bootstrapped_values.append(bt_df.corr(method=method).loc['exp_dG', column])

#                 bootstrapped_values = np.sort(bootstrapped_values)

#                 # Calculate confidence% confidence intervals
#                 lower = bootstrapped_values[int(tail * replicates)]
#                 upper = bootstrapped_values[int((1 - tail) * replicates)]

#                 if method == 'predictive_index':
#                     values[method][column] = (
#                         df[['exp_dG', column]].corr(method=predictive_index).loc['exp_dG', column],
#                         lower,
#                         upper,
#                     )
#                 else:
#                     values[method][column] = (
#                         df[['exp_dG', column]].corr(method=method).loc['exp_dG', column],
#                         lower,
#                         upper,
#                     )
#     out_df = pd.DataFrame(index=values[method].keys(), columns=values.keys())
#     for method in values:
#         for column in values[method]:
#             value, lower, upper = values[method][column]
#             # out_df.loc[column, method] = f"{value:.2f} [95%: {lower:.2f}, {upper:.2f}]"
#             out_df.loc[column, method] = f"{value:.2f} [{confidence}%: {lower:.2f}, {upper:.2f}]"
#     return out_df


def _bootstrap_correlation(df, method, column, replicates, tail):
    # Filter out rows with None or NaN in 'exp_dG' or the current simulation column
    filtered_df = df[['exp_dG', column]].dropna()

    # Perform bootstrapping
    bootstrapped_values = []
    for indexes in get_bootstrap_indexes(filtered_df, replicates=replicates):
        bt_df = filtered_df.iloc[indexes, :]

        bootstrapped_values.append(bt_df.corr(method=method).loc['exp_dG', column])

    bootstrapped_values = np.sort(bootstrapped_values)
    lower = bootstrapped_values[int(tail * replicates)]
    upper = bootstrapped_values[int((1 - tail) * replicates)]

    corr_value = filtered_df.corr(method=method).loc['exp_dG', column]

    return column, corr_value, lower, upper


def get_corr(df, replicates=1000, confidence=95):
    tail = (1 - confidence / 100) / 2
    values = dict()

    methods = ['pearson', 'kendall', 'spearman']
    with ProcessPoolExecutor() as executor:
        for method in methods:
            values[method] = dict()
            futures = [
                executor.submit(_bootstrap_correlation, df, method, column, replicates, tail)
                for column in df.columns
                if column.startswith('simulation_')
            ]
            for future in futures:
                column, corr_value, lower, upper = future.result()
                values[method][column] = (corr_value, lower, upper)

    # Prepare output DataFrame
    out_df = pd.DataFrame(index=values[method].keys(), columns=values.keys())
    for method in values:
        for column in values[method]:
            value, lower, upper = values[method][column]
            out_df.loc[column, method] = f"{value:.2f} [{confidence}%: {lower:.2f}, {upper:.2f}]"
    return out_df


def get_rmse(df, replicates=1000, confidence=95):
    columns_of_interest = [column for column in df.columns if column.startswith('simulation')]
    out_df = pd.DataFrame(index=columns_of_interest, columns=['rmse', 'mse', 'mue'])
    for column_of_interest in columns_of_interest:
        filtered_df = df[['exp_dG', column_of_interest]].dropna()
        intervals = confidence_intervals(
            functions=[rmse, mse, mue],
            actual_values=filtered_df['exp_dG'],
            predicted_values=filtered_df[column_of_interest],
            replicates=replicates,
            confidence=confidence,
        )
        rmse_value = (
            rmse(actual_values=filtered_df['exp_dG'], predicted_values=filtered_df[column_of_interest]),
            intervals[0][0],
            intervals[0][1],

        )
        rmse_value = [round(item, 2) for item in rmse_value]
        out_df.loc[column_of_interest, 'rmse'] = f"{rmse_value[0]} [{confidence}%: {rmse_value[1]}, {rmse_value[2]}]"

        mse_value = (
            mse(actual_values=filtered_df['exp_dG'], predicted_values=filtered_df[column_of_interest]),
            intervals[1][0],
            intervals[1][1],
        )
        mse_value = [round(item, 2) for item in mse_value]
        out_df.loc[column_of_interest, 'mse'] = f"{mse_value[0]} [{confidence}%: {mse_value[1]}, {mse_value[2]}]"

        mue_value = (
            mue(actual_values=filtered_df['exp_dG'], predicted_values=filtered_df[column_of_interest]),
            intervals[2][0],
            intervals[2][1],
        )
        mue_value = [round(item, 2) for item in mue_value]
        out_df.loc[column_of_interest, 'mue'] = f"{mue_value[0]} [{confidence}%: {mue_value[1]}, {mue_value[2]}]"

    return out_df


def predictive_index(data1, data2):
    """Implementation of predictive index from:
    https://pubs.acs.org/doi/epdf/10.1021/jm0100279
    The equations followed doi: 10.3389/fbinf.2022.885983, pg 15
    """
    def get_c_ij(delta_y, delta_x):
        if delta_x == 0:
            return 0
        elif delta_y / delta_x <= 0:
            return 1
        else:
            return 0
    assert len(data1) == len(data2)
    numerator = 0
    denominator = 0
    for i in range(len(data1)):
        for j in range(i+1, len(data2)):
            delta_y = data2[j] - data1[j]
            delta_x = data1[j] - data1[i]
            w_ij = abs(delta_y)
            c_ij = get_c_ij(delta_y, delta_x)

            numerator += w_ij*c_ij
            denominator += w_ij
    return numerator / denominator


def rmse(actual_values, predicted_values):
    actual_values = np.array(actual_values)
    predicted_values = np.array(predicted_values)
    return np.sqrt(np.mean((predicted_values - actual_values)**2))


def mue(actual_values, predicted_values):
    actual_values = np.array(actual_values)
    predicted_values = np.array(predicted_values)
    return np.mean(np.abs(predicted_values - actual_values))


def mse(actual_values, predicted_values):
    actual_values = np.array(actual_values)
    predicted_values = np.array(predicted_values)
    return np.mean(predicted_values - actual_values)


def get_bootstrap_indexes(sample: Iterable, replicates: int = 1000):
    # Perform bootstrapping
    bootstrap_indices = []
    for _ in range(replicates):
        # Generate a bootstrap sample by sampling with replacement
        bootstrap_indices.append(np.random.choice(len(sample), size=len(sample), replace=True))
    return bootstrap_indices


def confidence_intervals(functions: Union[Callable, Iterable[Callable]], actual_values: Iterable, predicted_values: Iterable,
                         replicates: int = 1000, confidence: float = 95) -> List[Tuple[float]]:
    """Calculate confidence intervals by bootstrapping with replacement

    Parameters
    ----------
    functions : Union[Callable, Iterable[Callable]]
        A list of functions to evaluate the bootstrapped data
    actual_values : Iterable
        Th reference values
    predicted_values : Iterable
        Predicted values
    replicates : int, optional
        Number of bootstrap replicates, by default 1000
    confidence : float, optional
        the percent of confidence, by default 95

    Returns
    -------
    List[Tuple[float]]
        A list with min and max for the interval in the order of functions
    """
    actual_values = np.array(actual_values)
    predicted_values = np.array(predicted_values)

    if not isinstance(functions, list):
        functions = [functions]

    values = [np.zeros(replicates) for _ in range(len(functions))]

    # Perform bootstrapping
    for i in range(replicates):
        # Generate a bootstrap sample by sampling with replacement
        bootstrap_indices = np.random.choice(len(actual_values), size=len(actual_values), replace=True)
        bootstrap_actual = actual_values[bootstrap_indices]
        bootstrap_predicted = predicted_values[bootstrap_indices]

        # Calculate function value the bootstrap sample
        for value, function in zip(values, functions):
            value[i] = function(bootstrap_actual, bootstrap_predicted)

    intervals = []
    for value in values:
        # Sort values
        sorted_values = np.sort(value)
        # Calculate confidence% confidence intervals
        tail = (1 - confidence / 100) / 2
        intervals.append(
            (
                sorted_values[int(tail * replicates)],
                sorted_values[int((1 - tail) * replicates)]
            )
        )
    return intervals
    return intervals

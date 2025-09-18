import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from numpy.polynomial.polynomial import polyfit
from utility import get_all_stats, get_numerical_values, fmt


def get_all_data(CALC_TYPE, FORCE_FIELD_ORDER, SYSTEM_ORDER, FILTER_THROMBIN_LIG_4):
    BindFlowData = pd.read_csv("../data/simulation/bindflow/gather/BindFlow.csv", index_col=0)
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
    if FILTER_THROMBIN_LIG_4:
        BindFlowData = BindFlowData[~((BindFlowData["source"] == "thrombin") & (BindFlowData["ligand"] == "lig_4"))]

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

    df_mse_removed = df_merge.copy()
    for ff in FORCE_FIELD_ORDER:
        for system_name in SYSTEM_ORDER:
            system_index = df_merge["source"] == system_name

            mse = (df_merge.loc[system_index, f"simulation_{ff}"] - df_merge.loc[system_index, "exp_dG"]).mean()
            df_mse_removed.loc[system_index, f"simulation_{ff}"] -= mse

            # uncertainty of the difference
            difference_uncertainty = np.sqrt(df_merge.loc[system_index, "exp_dG_error"]**2 + df_merge.loc[system_index, f"sem_{ff}"]**2)
            # uncertainty of the mean
            mse_uncertainty = np.sqrt(np.sum(difference_uncertainty**2)) / len(system_index)
            # uncertainty of the difference
            df_mse_removed.loc[system_index, f"sem_{ff}"] = np.sqrt(df_merge.loc[system_index, f"sem_{ff}"]**2 + mse_uncertainty**2)

    return df_merge, df_mse_removed


def make_publication_plot3(ax, df, simulation, system_name,
                           replicates=10_000, confidence=68,
                           min_value=-12, max_value=0,
                           main_fig=True,
                           remove_xticks_labels=False,
                           remove_yticks_labels=False,
                           return_stats=False
                           ):


    ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
    all_stats = get_all_stats(df, replicates=replicates, confidence=confidence)

    system_name = SYSTEM_NAME[system_name]
    sys_info_text = (
        r"~\\"
        fr"\textbf{{{system_name}}} \\[0.2em]"
        fr"$\mathrm{{N}} = {{{len(df)}}}$"
    )
    pearson = fmt(get_numerical_values(all_stats.loc[simulation, 'pearson']))
    kendall = fmt(get_numerical_values(all_stats.loc[simulation, 'kendall']))
    mue = fmt(get_numerical_values(all_stats.loc[simulation, 'mue']))
    rmse = fmt(get_numerical_values(all_stats.loc[simulation, 'rmse']))

    stats_text = (
        r"\begin{equation*}"
        r"\begin{aligned}[t]"
        fr"  \rho &= {pearson[0]}_{{{pearson[1]}}}^{{{pearson[2]}}} \\"
        fr"  \tau &= {kendall[0]}_{{{kendall[1]}}}^{{{kendall[2]}}} \\"
        fr"  \mathrm{{ocMUE}} &= {mue[0]}_{{{mue[1]}}}^{{{mue[2]}}} \\"
        fr"  \mathrm{{ocRMSE}} &= {rmse[0]}_{{{rmse[1]}}}^{{{rmse[2]}}}"
        r"\end{aligned}"
        r"\end{equation*}"
    )

    ax.text(
        0.02, 0.98, sys_info_text, transform=ax.transAxes, fontsize=35 if main_fig else 16,
        verticalalignment='top', horizontalalignment='left',
        )
    ax.text(
        0.98, 0.02, stats_text, transform=ax.transAxes, fontsize=31 if main_fig else 12,
        verticalalignment='bottom', horizontalalignment='right',
        )

    simulation_name = simulation.split('simulation_')[-1]
    # Plot the y = x line
    ax.plot([min_value, max_value], [min_value, max_value], '-', label='Y = X', color='black')
    # Fill the region between the error curves
    error = 2
    ax.fill_between([min_value, max_value],
                    [min_value - error, max_value - error],
                    [min_value + error, max_value + error],
                    alpha=0.3, label=r'$\pm$ 2 kcal/mol', color='gray', edgecolor='none')
    error = 1
    ax.fill_between([min_value, max_value],
                    [min_value - error, max_value - error],
                    [min_value + error, max_value + error],
                    alpha=0.3, label=r'$\pm$ 1 kcal/mol', color='gray', edgecolor='none')


    mask = df["exp_dG"].notna() & df[simulation].notna()
    if mask.sum() > 1:  # need at least 2 points to fit a line
        b, m = polyfit(df.loc[mask, "exp_dG"], df.loc[mask, simulation], 1)
        inter_range = np.array([min_value, max_value])
        ax.plot(inter_range, b + m * inter_range, '--', color='black')
    else:
        print(f"Not enough valid points for regression in {simulation}")

    # Plot the data

    # Scatter plot with color scale
    scatter = ax.scatter(
        df.exp_dG, df[simulation],
        c=(df.exp_dG - df[simulation]).abs(),
        cmap=plt.cm.coolwarm,
        norm=plt.Normalize(0, 5),
        s=100 if main_fig else 40, edgecolor='black', zorder=2)  # Higher zorder for dots

    # Add error bars behind the scatter points
    ax.errorbar(
        df.exp_dG, df[simulation],
        fmt='none', ecolor='gray',
        xerr=df.exp_dG_error, yerr=df[f"sem_{simulation_name}"], zorder=1)

    if not remove_yticks_labels:
        ax.set_ylabel(
            r"$\Delta G^{\mathrm{oc}}_{{\mathrm{{calc}}}} \, \mathrm{{[kcal/mol]}}$",
            fontsize=28 if main_fig else 16,
            ha='center',
        )

    if not remove_xticks_labels:
        ax.set_xlabel(
            r"$\Delta G_{\mathrm{exp}} \mathrm{[kcal/mol]}$",
            fontsize=28 if main_fig else 16)

    # Set equal scaling for x and y axes
    ax.set_aspect('equal', adjustable='box')

    ax.set_xlim([min_value, max_value])
    ax.set_ylim([min_value, max_value])

    # Define the tick values based on the common limits
    pad = 1
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
        start + 3 * step
    ]

    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    ax.tick_params(axis='x', labelsize=28 if main_fig else 15)
    ax.tick_params(axis='y', labelsize=28 if main_fig else 15)

    if remove_xticks_labels:
        ax.set_xticklabels([])  # Remove X-axis tick labels
    elif remove_yticks_labels:
        ax.set_yticklabels([])  # Remove Y-axis tick labels
    
    if return_stats:
        return scatter, (pearson, kendall, mue, rmse)
    else:
        return scatter


if __name__ == "__main__":
    # Warm-up hack to fix LaTeX rendering issue
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'  # Use Computer Modern fonts
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath,amsfonts}'
    plt.figure()
    plt.plot([0, 1], [0, 1])
    plt.xlabel(r"$\Delta G_{\mathrm{exp}}$")
    plt.ylabel(r"$\Delta G^{*}_{\mathrm{calc}}$")
    plt.close()

    replicates = 10_000
    confidence = 68

    FORCE_FIELD_ORDER = ['espaloma-0.3.1', 'gaff-2.11', 'openff-2.0.0']
    SYSTEM_ORDER = ["p38", "A2A", "ptp1b", "tyk2", "thrombin", "mcl1", "CyclophilinD"]
    SYSTEM_NAME = {
        "All": "All",
        "p38": "P38",
        "A2A": "A2A",
        "ptp1b": "PTP1B",
        "tyk2": "TYK2",
        "thrombin": "Thrombin",
        "mcl1": "MCL1",
        "CyclophilinD": "CyclophilinD",
        "SAMPL6-OA": "SAMPL6-OA"
    }
    collapsed_metrics = []
    for CALC_TYPE in ["mbar", "dh-gb"]:
        for FILTER_THROMBIN_LIG_4 in [True, False]:
            # Plotting
            for SELECTED_FORCE_FIELD in FORCE_FIELD_ORDER:
                print(CALC_TYPE, FILTER_THROMBIN_LIG_4, SELECTED_FORCE_FIELD)
                local_system_order = SYSTEM_ORDER.copy()
                if SELECTED_FORCE_FIELD == 'espaloma-0.3.1':
                    local_system_order.append("SAMPL6-OA")

                _, df_merge_mse_removed = get_all_data(
                    CALC_TYPE=CALC_TYPE,
                    FORCE_FIELD_ORDER=FORCE_FIELD_ORDER,
                    SYSTEM_ORDER=local_system_order,
                    FILTER_THROMBIN_LIG_4=FILTER_THROMBIN_LIG_4
                    )

                if SELECTED_FORCE_FIELD == 'espaloma-0.3.1':
                    fig = plt.figure(figsize=(16, 18))
                    gs = GridSpec(5, 4, figure=fig, hspace=0, wspace=0.1)
                    ax1 = fig.add_subplot(gs[1:4, 0:3])

                    axb1 = fig.add_subplot(gs[4, 0])
                    axb2 = fig.add_subplot(gs[4, 1])
                    axb3 = fig.add_subplot(gs[4, 2])
                    axb4 = fig.add_subplot(gs[4, 3])

                    axr1 = fig.add_subplot(gs[0, 3])
                    axr2 = fig.add_subplot(gs[1, 3])
                    axr3 = fig.add_subplot(gs[2, 3])
                    axr4 = fig.add_subplot(gs[3, 3])
                else:
                    fig = plt.figure(figsize=(16, 14))
                    gs = GridSpec(4, 4, figure=fig, hspace=0, wspace=0.1)

                    ax1 = fig.add_subplot(gs[0:3, 0:3])
                    axb1 = fig.add_subplot(gs[3, 0])
                    axb2 = fig.add_subplot(gs[3, 1])
                    axb3 = fig.add_subplot(gs[3, 2])
                    axb4 = fig.add_subplot(gs[3, 3])

                    axr1 = fig.add_subplot(gs[0, 3])
                    axr2 = fig.add_subplot(gs[1, 3])
                    axr3 = fig.add_subplot(gs[2, 3])

                box = ax1.get_position()
                ax1.set_position([box.x0, box.y0 + 0.05, box.width * 0.9, box.height * 0.9])
                # last row
                if SELECTED_FORCE_FIELD == 'espaloma-0.3.1':
                    axis_individual = [axb1, axb2, axb3, axb4, axr4, axr3, axr2, axr1]
                else:
                    axis_individual = [axb1, axb2, axb3, axb4, axr3, axr2, axr1]

                padding = 1
                min_value = np.min(df_merge_mse_removed[[f"simulation_{SELECTED_FORCE_FIELD}", "exp_dG"]]) - padding
                max_value = np.max(df_merge_mse_removed[[f"simulation_{SELECTED_FORCE_FIELD}", "exp_dG"]]) + padding

                for ax, system in zip(axis_individual, local_system_order):
                    if system in ["p38", "A2A", "ptp1b", "tyk2"]:
                        remove_xticks_labels = False
                        remove_yticks_labels = True
                        if system == "p38":
                            remove_yticks_labels = False
                    elif system in ["SAMPL6-OA", "CyclophiliD", "mcl1", "thrombin"]:
                        remove_xticks_labels = True
                        remove_yticks_labels = False
                    _ = make_publication_plot3(
                        ax, df=df_merge_mse_removed[df_merge_mse_removed["source"] == system], simulation=f"simulation_{SELECTED_FORCE_FIELD}",
                        system_name=system, replicates=replicates, confidence=confidence, min_value=min_value,
                        max_value=max_value, main_fig=False,
                        remove_xticks_labels=remove_xticks_labels,
                        remove_yticks_labels=remove_yticks_labels
                    )

                scatter, stats = make_publication_plot3(
                    ax1, df=df_merge_mse_removed, simulation=f"simulation_{SELECTED_FORCE_FIELD}",
                    system_name="All", replicates=replicates, confidence=confidence,
                    min_value=min_value,
                    max_value=max_value,
                    return_stats=True
                    )

                collapsed_metrics.append([CALC_TYPE, not FILTER_THROMBIN_LIG_4, SELECTED_FORCE_FIELD, *stats])

                # Add colorbar
                if SELECTED_FORCE_FIELD == "espaloma-0.3.1":
                    shrink = 0.63
                else:
                    shrink = 0.8
                cbar = fig.colorbar(scatter, ax=[ax1] + axis_individual, location="right", shrink=shrink, pad=0.01, aspect=30)  # Attach the colorbar to the current axis
                cbar.set_label(r'$\lvert \Delta G^{\mathrm{oc}}_{\mathrm{calc}} - \Delta G_{\mathrm{exp}} \rvert$ [kcal/mol]', fontsize=28)
                cbar.ax.tick_params(labelsize=28)
                # fig.subplots_adjust(wspace=0.1, hspace=0.1)  # Adjust the value as needed
                if FILTER_THROMBIN_LIG_4:
                    fig.savefig(f"summary/ff/{CALC_TYPE}_{SELECTED_FORCE_FIELD}_thrombin_lig_4_removed.pdf",
                                bbox_inches="tight",
                                pad_inches=0.05,
                                transparent=False)
                else:
                    fig.savefig(f"summary/ff/{CALC_TYPE}_{SELECTED_FORCE_FIELD}.pdf",
                                bbox_inches="tight",
                                pad_inches=0.05,
                                transparent=False)

    # pd.DataFrame(
    #     columns=["calc_type", "with_outlier", "ff", "pearson", "kendall", "mue", "rmse"],
    #     data=collapsed_metrics
    # ).to_csv("summary/ff/collapsed_metrics.csv")

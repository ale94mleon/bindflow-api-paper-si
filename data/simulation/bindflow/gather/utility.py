from typing import Callable, Iterable, List, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython.display import display
from numpy.polynomial.polynomial import polyfit
from plbenchmark import edges, ligands, targets
from concurrent.futures import ProcessPoolExecutor




def get_all_stats(df, replicates=1000, confidence=95):
    all_stats = pd.merge(
        get_corr(df, replicates=replicates, confidence=confidence),
        get_rmse(df, replicates=replicates, confidence=confidence),
        left_index=True,
        right_index=True)
    return all_stats


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

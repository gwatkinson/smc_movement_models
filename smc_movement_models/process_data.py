"""Module pour nettoyer les données.

Usage:
    python smc_movement_models/process_data.py
"""

import click
import pandas as pd

RAW_DATA_PATH = "./data/4121_filtervelocity_PostThesis.csv"
DAY = "18/01/2008"
DAY_PATH = "./data/18_01_2008.csv"
OUTPUT_TOTAL_PATH = "./data/clean_data.csv"
SUMMARY_PATH = "./data/summary.csv"


@click.command()
@click.option(
    "-p", "--path", default=RAW_DATA_PATH, help="Path to raw data", prompt="Path to original data"
)
@click.option(
    "-d",
    "--day",
    default=DAY,
    help="Day to keep",
    prompt="Which day to keep (from 14/01/2008 to 28/01/2008) ?",
)
@click.option(
    "-o",
    "--day-path",
    default=DAY_PATH,
    help="Path to output day data",
    prompt="Name of output day file",
)
@click.option(
    "-t",
    "--output-total",
    default=OUTPUT_TOTAL_PATH,
    help="Path to output total data",
    prompt="Name of output clean file",
)
@click.option(
    "-s",
    "--summary-path",
    default=SUMMARY_PATH,
    help="Path to output summary data",
    prompt="Name of output summary file",
)
def clean_raw_data(path, day, output_total, summary_path, day_path):
    # Raw data
    click.echo(f"Reading raw data from: {path}")
    df = pd.read_csv(path)
    df["Dtime"] = pd.to_datetime(df["Dtime"])
    df["Depth_Diff"] = df["Calib_Depth"] - df["Orig_Depth"]
    df["Depth"] = -df["Velocity"].cumsum() * 5
    df["Cumulative_Sum_Diff"] = -df["Depth"] - df["Calib_Depth"]

    click.echo(f"Saving total data to: {output_total}")
    df.to_csv(output_total, index=False)

    # Day data
    mask_day = df["Date"] == day
    day_df = df.loc[mask_day, ["Dtime", "Velocity", "Depth"]]
    click.echo(f"Saving day data to: {day_path}")
    day_df.to_csv(day_path, index=False)

    # Summary
    nan_sum = lambda x: x.isnull().sum()
    summary = df.groupby("Date").agg(
        Min_Time=pd.NamedAgg(column="Dtime", aggfunc="min"),
        Max_Time=pd.NamedAgg(column="Dtime", aggfunc="max"),
        Calib_Depth_NaNs=pd.NamedAgg(column="Calib_Depth", aggfunc=nan_sum),
        Orig_Depth_NaNs=pd.NamedAgg(column="Orig_Depth", aggfunc=nan_sum),
        Velocity_NaNs=pd.NamedAgg(column="Velocity", aggfunc=nan_sum),
        Observations=pd.NamedAgg(column="Velocity", aggfunc="count"),
    )
    click.echo(f"Saving summary to: {summary_path}")
    summary.to_csv(summary_path, index=False)


if __name__ == "__main__":
    clean_raw_data()

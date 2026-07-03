import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

matplotlib.use("Agg")
pd.options.mode.chained_assignment = None


def plot_prc(accuracy_file, output_file, plot_configs):

    df = pd.read_csv(accuracy_file, sep="\t").dropna()

    if "threshold" in df.columns:
        df.rename(columns={"threshold": "cutoff"}, inplace=True)

    # columns = ['demography', 'scenario', 'replicate', 'cutoff', 'sample']
    columns = ["demography", "scenario", "cutoff", "sample"]

    columns = [x for x in columns if x in df.columns]

    matched_model = ""

    for m in plot_configs["models_for_plotting"]:
        if m in accuracy_file:
            matched_model = m
            break

    recall = df["recall"]
    precision = df["precision"]

    fig, axs = plt.subplots(
        nrows=2, ncols=3, constrained_layout=True, figsize=(7.5, 4), dpi=350
    )
    gridspec = axs[0, 0].get_subplotspec().get_gridspec()
    for a in axs[:, 2]:
        a.remove()
    j = 0

    # average
    # Group by relevant columns
    df_mean = df.groupby(columns, as_index=False).mean()

    def f1_iso(recall, f1=0.5):
        """Return precision for a given F1 score and recall."""
        return (f1 * recall) / (2 * recall - f1 + 1e-8)

    plt.figure(figsize=(5, 4))

    # Plot PR curve
    for scenario in df_mean["scenario"].unique():
        subset = df_mean[df_mean["scenario"] == scenario]
        plt.plot(subset["recall"], subset["precision"], marker="o", label=scenario)

    # F1 iso-curves
    recall_vals = np.linspace(0.01, 1, 100)
    f_scores = np.linspace(20, 80, num=4)
    for f1 in f_scores:
        plt.plot(
            recall_vals,
            f1_iso(recall_vals, f1=f1),
            color="gray",
            linestyle=":",
            alpha=0.5,
        )
        plt.text(
            0.95, f1_iso(0.95, f1), f"F1={f1}", fontsize=8, va="bottom", ha="right"
        )

    for f_score in f_scores:
        x = np.linspace(1, 100)
        y = f_score * x / (2 * x - f_score)
        (l,) = axs[0, j].plot(
            x[y >= 0], y[y >= 0], color="black", alpha=0.4, linestyle="dotted", zorder=1
        )
        axs[0, j].annotate(
            "F1={0:0.0f}%".format(f_score), xy=(101, y[45] + 2), fontsize=8
        )

    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title(str(matched_model) + ": Precision–Recall Curve")
    plt.legend()
    plt.tight_layout()
    plt.show()

    plt.savefig(output_file, bbox_inches="tight")


def plot_confusion_matrices(
    file_path, output_pdf, thresholds_to_plot=None, best_f1_only=False
):
    """
        Plots confusion matrices for scenarios in a tab-separated file and saves them to a PDF.

        Positive class = 'introgressed'
        Negative class = 'not introgressed'
    :
                  Predicted
                  not introgressed   introgressed
        Actual
        not introgressed      TN           FP
        introgressed          FN           TP

        Parameters:
        - file_path: str, path to the .txt file
        - output_pdf: str, path of the PDF file to save plots
        - thresholds_to_plot: list of ints or floats, optional thresholds to include
        - best_f1_only: bool, if True, plot only the threshold(s) with the highest F1 per scenario
    """
    from matplotlib.backends.backend_pdf import PdfPages

    df = pd.read_csv(file_path, sep="\t")

    if "cutoff" in df.columns and "threshold" not in df.columns:
        df = df.rename(columns={"cutoff": "threshold"})

    # filter by thresholds
    if thresholds_to_plot is not None:
        df = df[df["threshold"].isin(thresholds_to_plot)]

    # Group by scenario
    grouped = df.groupby("scenario")

    # Open PDF
    with PdfPages(output_pdf) as pdf:
        for scenario, group in grouped:
            if best_f1_only:
                max_f1 = group["F1"].max()
                group = group[group["F1"] == max_f1]

            for _, row in group.iterrows():
                tn, fp, fn, tp = row["TN"], row["FP"], row["FN"], row["TP"]
                threshold = row["threshold"]
                title = f'{scenario} | Threshold={threshold} | F1={row["F1"]:.2f}'

                cm = [[tn, fp], [fn, tp]]

                plt.figure(figsize=(6, 5))
                sns.heatmap(
                    cm,
                    annot=True,
                    fmt="g",
                    cmap="Blues",
                    xticklabels=["not introgressed", "introgressed"],
                    yticklabels=["not introgressed", "introgressed"],
                )
                plt.xlabel("Predicted")
                plt.ylabel("Actual")
                plt.title(title)

                pdf.savefig()
                plt.close()

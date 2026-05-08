#!/bin/env python

import argparse
import re
import textwrap

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
import pandas as pd
import seaborn as sns


def pdf_output(path):
    if not path.lower().endswith(".pdf"):
        raise argparse.ArgumentTypeError('--output must end with ".pdf"')
    return path


def optional_float(value):
    if value is None:
        return None

    value = str(value).lower()
    if value == "none":
        return None

    return float(value)


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument(
    "--input",
    "-i",
    type=str,
    required=True,
    help="all_pairwise.tsv from sample-id-check",
)

parser.add_argument(
    "--output",
    "-o",
    type=pdf_output,
    required=True,
    help='output PDF name. It should end with ".pdf".',
)

parser.add_argument(
    "--preview-png",
    type=str,
    default=None,
    help="Optional PNG preview output path",
)

parser.add_argument(
    "--annot-max-samples",
    type=int,
    default=25,
    help="Show cell annotations only when sample count is <= this value",
)

parser.add_argument(
    "--max-figsize",
    type=optional_float,
    default=22,
    help='Maximum figure size in inches. Use "none" for no limit.',
)

parser.add_argument(
    "--cell-size",
    type=float,
    default=1.1,
    help="Figure size per sample. Increase this to make heatmap larger.",
)

parser.add_argument(
    "--max-labels",
    type=int,
    default=80,
    help="Maximum number of x/y labels to show before thinning labels",
)

parser.add_argument(
    "--title",
    type=str,
    default="Data Distance Heatmap",
    help="Heatmap title",
)

args = parser.parse_args()


def build_distance_matrix(df):
    samples = sorted(set(df["sample1"]).union(df["sample2"]))

    distance_matrix = pd.DataFrame(
        np.inf,
        index=samples,
        columns=samples,
    )

    for _, row in df.iterrows():
        sample1 = row["sample1"]
        sample2 = row["sample2"]
        score = row["score"]

        distance_matrix.loc[sample1, sample2] = score
        distance_matrix.loc[sample2, sample1] = score

    np.fill_diagonal(distance_matrix.values, 0)

    return distance_matrix, samples


def build_colormap(data_max):
    vmax = max(4.0, data_max)

    cool_blue = "#3b4cc0"
    light_blue = "#82a6d3"
    middle_color = "#f7f7f7"
    warm_orange = "#f07c52"
    warm_red = "#b40426"

    colors = [
        (0.0, cool_blue),
        (0.5, light_blue),
        (1.0, middle_color),
        (min(max(data_max, 2.0), 4.0), warm_orange),
        (max(data_max, 4.0), warm_red),
    ]

    color_stops = [
        (min(v / vmax, 1.0), color)
        for v, color in sorted(colors, key=lambda x: x[0])
    ]

    custom_cmap = LinearSegmentedColormap.from_list(
        "data_distance_map",
        color_stops,
    )

    norm = Normalize(vmin=0.0, vmax=vmax)

    return custom_cmap, norm


def clean_label(label):
    label = str(label)
    label = label.replace(".count", "")
    label = re.sub(r"\s+", " ", label)
    return label


def wrap_label(label, max_line_len=28, max_lines=3):
    label = clean_label(label)

    if len(label) <= max_line_len:
        return label

    # Prefer line breaks at natural separators.
    label_for_wrap = re.sub(r"([/_-])", r"\1 ", label)

    wrapped_lines = textwrap.wrap(
        label_for_wrap,
        width=max_line_len,
        break_long_words=False,
        break_on_hyphens=True,
    )

    cleaned_lines = [
        line.replace(" ", "")
        for line in wrapped_lines
    ]

    if len(cleaned_lines) <= max_lines:
        return "\n".join(cleaned_lines)

    kept_lines = cleaned_lines[:max_lines]
    kept_lines[-1] = kept_lines[-1][: max_line_len - 3] + "..."

    return "\n".join(kept_lines)


def get_label_wrap_width(n_samples):
    if n_samples <= 20:
        return 36
    if n_samples <= 50:
        return 30
    if n_samples <= 100:
        return 24
    return 18


def get_label_max_lines(n_samples):
    if n_samples <= 30:
        return 4
    if n_samples <= 80:
        return 3
    return 2


def get_tick_step(n_samples):
    if n_samples <= args.max_labels:
        return 1

    return int(np.ceil(n_samples / args.max_labels))


def get_font_sizes(n_samples):
    tick_size = max(
        5,
        min(13, int(260 / max(n_samples, 1))),
    )

    annot_size = max(
        5,
        min(12, int(220 / max(n_samples, 1))),
    )

    title_size = max(18, tick_size + 11)

    return tick_size, annot_size, title_size


def get_figure_size(samples, n_samples):
    cleaned_labels = [clean_label(sample) for sample in samples]
    max_label_len = max(len(label) for label in cleaned_labels)

    heatmap_size = max(12, n_samples * args.cell_size)

    # Label space is added outside of heatmap size.
    # This keeps plot large even when labels are long.
    label_space = min(max(5, max_label_len * 0.10), 18)

    fig_width = heatmap_size + label_space + 4
    fig_height = heatmap_size + label_space + 4

    if args.max_figsize is not None:
        fig_width = min(fig_width, args.max_figsize)
        fig_height = min(fig_height, args.max_figsize)

    return fig_width, fig_height


def should_show_annotation(n_samples):
    return n_samples <= args.annot_max_samples


def set_axis_labels(ax, samples, n_samples, tick_size):
    wrap_width = get_label_wrap_width(n_samples)
    max_lines = get_label_max_lines(n_samples)

    labels = [
        wrap_label(
            sample,
            max_line_len=wrap_width,
            max_lines=max_lines,
        )
        for sample in samples
    ]

    step = get_tick_step(n_samples)

    tick_positions = np.arange(0, n_samples, step) + 0.5
    tick_labels = labels[::step]

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(
        tick_labels,
        rotation=90,
        fontsize=tick_size,
        ha="center",
        va="top",
    )

    ax.set_yticks(tick_positions)
    ax.set_yticklabels(
        tick_labels,
        rotation=0,
        fontsize=tick_size,
        ha="right",
        va="center",
    )


def get_margins(fig_width, fig_height, samples, n_samples):
    max_label_len = max(len(clean_label(sample)) for sample in samples)

    label_pressure = min(max_label_len / 120, 0.18)

    left_margin = max(0.10, min(0.28, 4.0 / fig_width + label_pressure))
    bottom_margin = max(0.12, min(0.30, 4.5 / fig_height + label_pressure))

    # For very large plots, margins should not consume too much heatmap space.
    if n_samples > 80:
        left_margin = min(left_margin, 0.18)
        bottom_margin = min(bottom_margin, 0.20)

    return {
        "left": left_margin,
        "bottom": bottom_margin,
        "right": 0.94,
        "top": 0.91,
    }


def check_plot_layout(samples, n_samples, fig_width, fig_height, show_annot):
    max_label_len = max(len(clean_label(sample)) for sample in samples)
    warnings = []

    if fig_width < 12:
        warnings.append("Figure width is small.")

    if fig_height < 12:
        warnings.append("Figure height is small.")

    if n_samples * args.cell_size < 10:
        warnings.append("Heatmap area may be small. Consider larger --cell-size.")

    if show_annot and n_samples > 25:
        warnings.append("Annotations may be crowded.")

    if max_label_len > 80:
        warnings.append("Some labels are very long and will be wrapped/truncated.")

    if n_samples > 100:
        warnings.append("Large matrix detected. Labels will be thinned.")

    if warnings:
        print("\nPlot layout warnings:")
        for warning in warnings:
            print(f"  - {warning}")
    else:
        print("\nPlot layout check passed.")


def plot_heatmap(distance_matrix, samples, output):
    n_samples = len(samples)

    fig_width, fig_height = get_figure_size(samples, n_samples)
    tick_size, annot_size, title_size = get_font_sizes(n_samples)
    show_annot = should_show_annotation(n_samples)

    data_values = distance_matrix.replace(np.inf, np.nan).values
    data_max = np.nanmax(data_values)

    custom_cmap, norm = build_colormap(data_max)

    fig, ax = plt.subplots(
        figsize=(fig_width, fig_height),
        constrained_layout=False,
    )

    sns.heatmap(
        distance_matrix,
        annot=show_annot,
        fmt=".2f",
        cmap=custom_cmap,
        annot_kws={"size": annot_size},
        cbar=True,
        cbar_kws={
            "shrink": 0.75,
            "pad": 0.02,
        },
        norm=norm,
        square=False,
        linewidths=0.15 if n_samples <= 60 else 0,
        linecolor="white",
        ax=ax,
    )

    # Keep PDF size and rendering time manageable.
    ax.collections[0].set_rasterized(True)

    ax.set_title(
        args.title,
        fontsize=title_size,
        pad=20,
    )

    ax.set_xlabel("")
    ax.set_ylabel("")

    set_axis_labels(
        ax=ax,
        samples=samples,
        n_samples=n_samples,
        tick_size=tick_size,
    )

    margins = get_margins(
        fig_width=fig_width,
        fig_height=fig_height,
        samples=samples,
        n_samples=n_samples,
    )

    fig.subplots_adjust(**margins)

    check_plot_layout(
        samples=samples,
        n_samples=n_samples,
        fig_width=fig_width,
        fig_height=fig_height,
        show_annot=show_annot,
    )

    fig.savefig(
        output,
        format="pdf",
        dpi=300,
        bbox_inches="tight",
    )

    if args.preview_png is not None:
        fig.savefig(
            args.preview_png,
            format="png",
            dpi=200,
        )

    plt.close(fig)


def main():
    df = pd.read_csv(
        args.input,
        sep="\t",
        header=0,
        usecols=["sample1", "sample2", "score"],
    )

    distance_matrix, samples = build_distance_matrix(df)

    plot_heatmap(
        distance_matrix=distance_matrix,
        samples=samples,
        output=args.output,
    )


if __name__ == "__main__":
    main()